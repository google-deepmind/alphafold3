# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for running the MSA and template tools for the AlphaFold model."""

from concurrent import futures
import dataclasses
import datetime
import functools
import logging
import time

from alphafold3.common import folding_input
from alphafold3.constants import mmcif_names
from alphafold3.data import msa
from alphafold3.data import msa_config
from alphafold3.data import structure_stores
from alphafold3.data import templates as templates_lib


# Cache to avoid re-running template search for the same sequence in homomers.
@functools.cache
def _get_protein_templates(
    sequence: str,
    input_msa_a3m: str,
    run_template_search: bool,
    templates_config: msa_config.TemplatesConfig,
    pdb_database_path: str,
) -> templates_lib.Templates:
  """Searches for templates for a single protein chain."""
  if run_template_search:
    templates_start_time = time.time()
    logging.info('Getting protein templates for sequence %s', sequence)
    protein_templates = templates_lib.Templates.from_seq_and_a3m(
        query_sequence=sequence,
        msa_a3m=input_msa_a3m,
        max_template_date=templates_config.filter_config.max_template_date,
        database_path=templates_config.template_tool_config.database_path,
        hmmsearch_config=templates_config.template_tool_config.hmmsearch_config,
        max_a3m_query_sequences=None,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        structure_store=structure_stores.StructureStore(pdb_database_path),
        filter_config=templates_config.filter_config,
    )
    logging.info(
        'Getting %d protein templates took %.2f seconds for sequence %s',
        protein_templates.num_hits,
        time.time() - templates_start_time,
        sequence,
    )
  else:
    logging.info('Skipping template search for sequence %s', sequence)
    protein_templates = templates_lib.Templates(
        query_sequence=sequence,
        hits=[],
        max_template_date=templates_config.filter_config.max_template_date,
        structure_store=structure_stores.StructureStore(pdb_database_path),
    )
  return protein_templates


# Cache to avoid re-running the MSA tools for the same sequence in homomers.
@functools.cache
def _get_protein_msa_and_templates(
    sequence: str,
    run_template_search: bool,
    uniref90_msa_config: msa_config.RunConfig,
    mgnify_msa_config: msa_config.RunConfig,
    small_bfd_msa_config: msa_config.RunConfig,
    uniprot_msa_config: msa_config.RunConfig,
    templates_config: msa_config.TemplatesConfig,
    pdb_database_path: str,
) -> tuple[msa.Msa, msa.Msa, templates_lib.Templates]:
  """Processes a single protein chain."""
  logging.info('Getting protein MSAs for sequence %s', sequence)
  msa_start_time = time.time()
  # Run various MSA tools in parallel. Use a ThreadPoolExecutor because
  # they're not blocked by the GIL, as they're sub-shelled out.
  with futures.ThreadPoolExecutor(max_workers=4) as executor:
    uniref90_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=uniref90_msa_config,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    mgnify_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=mgnify_msa_config,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    small_bfd_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=small_bfd_msa_config,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    uniprot_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=uniprot_msa_config,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
  uniref90_msa = uniref90_msa_future.result()
  mgnify_msa = mgnify_msa_future.result()
  small_bfd_msa = small_bfd_msa_future.result()
  uniprot_msa = uniprot_msa_future.result()
  logging.info(
      'Getting protein MSAs took %.2f seconds for sequence %s',
      time.time() - msa_start_time,
      sequence,
  )

  logging.info('Deduplicating MSAs for sequence %s', sequence)
  msa_dedupe_start_time = time.time()
  with futures.ThreadPoolExecutor() as executor:
    unpaired_protein_msa_future = executor.submit(
        msa.Msa.from_multiple_msas,
        msas=[uniref90_msa, small_bfd_msa, mgnify_msa],
        deduplicate=True,
    )
    paired_protein_msa_future = executor.submit(
        msa.Msa.from_multiple_msas, msas=[uniprot_msa], deduplicate=False
    )
  unpaired_protein_msa = unpaired_protein_msa_future.result()
  paired_protein_msa = paired_protein_msa_future.result()
  logging.info(
      'Deduplicating MSAs took %.2f seconds for sequence %s, found %d unpaired'
      ' sequences, %d paired sequences',
      time.time() - msa_dedupe_start_time,
      sequence,
      unpaired_protein_msa.depth,
      paired_protein_msa.depth,
  )

  protein_templates = _get_protein_templates(
      sequence=sequence,
      input_msa_a3m=unpaired_protein_msa.to_a3m(),
      run_template_search=run_template_search,
      templates_config=templates_config,
      pdb_database_path=pdb_database_path,
  )

  return unpaired_protein_msa, paired_protein_msa, protein_templates


# Cache to avoid re-running the Nhmmer for the same sequence in homomers.
@functools.cache
def _get_rna_msa(
    sequence: str,
    nt_rna_msa_config: msa_config.NhmmerConfig,
    rfam_msa_config: msa_config.NhmmerConfig,
    rnacentral_msa_config: msa_config.NhmmerConfig,
) -> msa.Msa:
  """Processes a single RNA chain."""
  logging.info('Getting RNA MSAs for sequence %s', sequence)
  rna_msa_start_time = time.time()
  # Run various MSA tools in parallel. Use a ThreadPoolExecutor because
  # they're not blocked by the GIL, as they're sub-shelled out.
  with futures.ThreadPoolExecutor() as executor:
    nt_rna_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=nt_rna_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
    rfam_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=rfam_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
    rnacentral_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=rnacentral_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
  nt_rna_msa = nt_rna_msa_future.result()
  rfam_msa = rfam_msa_future.result()
  rnacentral_msa = rnacentral_msa_future.result()
  rna_msa = msa.Msa.from_multiple_msas(
      msas=[rfam_msa, rnacentral_msa, nt_rna_msa],
      deduplicate=True,
  )
  logging.info(
      'Getting RNA MSAs took %.2f seconds for sequence %s, found %d unpaired'
      ' sequences',
      time.time() - rna_msa_start_time,
      sequence,
      rna_msa.depth,
  )
  return rna_msa


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class DataPipelineConfig:
  """The configuration for the data pipeline.

  Attributes:
    jackhmmer_binary_path: Jackhmmer binary path, used for protein MSA search.
    nhmmer_binary_path: Nhmmer binary path, used for RNA MSA search.
    hmmalign_binary_path: Hmmalign binary path, used to align hits to the query
      profile.
    hmmsearch_binary_path: Hmmsearch binary path, used for template search.
    hmmbuild_binary_path: Hmmbuild binary path, used to build HMM profile from
      raw MSA in template search.
    small_bfd_database_path: Small BFD database path, used for protein MSA
      search.
    mgnify_database_path: Mgnify database path, used for protein MSA search.
    uniprot_cluster_annot_database_path: Uniprot database path, used for protein
      paired MSA search.
    uniref90_database_path: UniRef90 database path, used for MSA search, and the
      MSA obtained by searching it is used to construct the profile for template
      search.
    ntrna_database_path: NT-RNA database path, used for RNA MSA search.
    rfam_database_path: Rfam database path, used for RNA MSA search.
    rna_central_database_path: RNAcentral database path, used for RNA MSA
      search.
    seqres_database_path: PDB sequence database path, used for template search.
    pdb_database_path: PDB database directory with mmCIF files path, used for
      template search.
    jackhmmer_n_cpu: Number of CPUs to use for Jackhmmer.
    nhmmer_n_cpu: Number of CPUs to use for Nhmmer.
    max_template_date: The latest date of templates to use.
  """

  # Binary paths.
  jackhmmer_binary_path: str
  nhmmer_binary_path: str
  hmmalign_binary_path: str
  hmmsearch_binary_path: str
  hmmbuild_binary_path: str

  # Jackhmmer databases.
  small_bfd_database_path: str
  mgnify_database_path: str
  uniprot_cluster_annot_database_path: str
  uniref90_database_path: str
  # Nhmmer databases.
  ntrna_database_path: str
  rfam_database_path: str
  rna_central_database_path: str
  # Template search databases.
  seqres_database_path: str
  pdb_database_path: str

  # Optional configuration for MSA tools.
  jackhmmer_n_cpu: int = 8
  nhmmer_n_cpu: int = 8

  max_template_date: datetime.date


class DataPipeline:
  """Runs the alignment tools and assembles the input features."""

  def __init__(self, data_pipeline_config: DataPipelineConfig):
    """Initializes the data pipeline with default configurations."""
    self._uniref90_msa_config = msa_config.RunConfig(
        config=msa_config.JackhmmerConfig(
            binary_path=data_pipeline_config.jackhmmer_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='uniref90',
                path=data_pipeline_config.uniref90_database_path,
            ),
            n_cpu=data_pipeline_config.jackhmmer_n_cpu,
            n_iter=1,
            e_value=1e-4,
            z_value=None,
            max_sequences=10_000,
        ),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        crop_size=None,
    )
    self._mgnify_msa_config = msa_config.RunConfig(
        config=msa_config.JackhmmerConfig(
            binary_path=data_pipeline_config.jackhmmer_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='mgnify',
                path=data_pipeline_config.mgnify_database_path,
            ),
            n_cpu=data_pipeline_config.jackhmmer_n_cpu,
            n_iter=1,
            e_value=1e-4,
            z_value=None,
            max_sequences=5_000,
        ),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        crop_size=None,
    )
    self._small_bfd_msa_config = msa_config.RunConfig(
        config=msa_config.JackhmmerConfig(
            binary_path=data_pipeline_config.jackhmmer_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='small_bfd',
                path=data_pipeline_config.small_bfd_database_path,
            ),
            n_cpu=data_pipeline_config.jackhmmer_n_cpu,
            n_iter=1,
            e_value=1e-4,
            # Set z_value=138_515_945 to match the z_value used in the paper.
            # In practice, this has minimal impact on predicted structures.
            z_value=None,
            max_sequences=5_000,
        ),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        crop_size=None,
    )
    self._uniprot_msa_config = msa_config.RunConfig(
        config=msa_config.JackhmmerConfig(
            binary_path=data_pipeline_config.jackhmmer_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='uniprot_cluster_annot',
                path=data_pipeline_config.uniprot_cluster_annot_database_path,
            ),
            n_cpu=data_pipeline_config.jackhmmer_n_cpu,
            n_iter=1,
            e_value=1e-4,
            z_value=None,
            max_sequences=50_000,
        ),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        crop_size=None,
    )
    self._nt_rna_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='nt_rna',
                path=data_pipeline_config.ntrna_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10_000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )
    self._rfam_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='rfam_rna',
                path=data_pipeline_config.rfam_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10_000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )
    self._rnacentral_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='rna_central_rna',
                path=data_pipeline_config.rna_central_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10_000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )

    self._templates_config = msa_config.TemplatesConfig(
        template_tool_config=msa_config.TemplateToolConfig(
            database_path=data_pipeline_config.seqres_database_path,
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
            hmmsearch_config=msa_config.HmmsearchConfig(
                hmmsearch_binary_path=data_pipeline_config.hmmsearch_binary_path,
                hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
                filter_f1=0.1,
                filter_f2=0.1,
                filter_f3=0.1,
                e_value=100,
                inc_e=100,
                dom_e=100,
                incdom_e=100,
                alphabet='amino',
            ),
        ),
        filter_config=msa_config.TemplateFilterConfig(
            max_subsequence_ratio=0.95,
            min_align_ratio=0.1,
            min_hit_length=10,
            deduplicate_sequences=True,
            max_hits=4,
            max_template_date=data_pipeline_config.max_template_date,
        ),
    )
    self._pdb_database_path = data_pipeline_config.pdb_database_path

  def process_protein_chain(
      self, chain: folding_input.ProteinChain
  ) -> folding_input.ProteinChain:
    """Processes a single protein chain."""
    has_unpaired_msa = chain.unpaired_msa is not None
    has_paired_msa = chain.paired_msa is not None
    has_templates = chain.templates is not None

    if not has_unpaired_msa and not has_paired_msa and not chain.templates:
      # MSA None - search. Templates either [] - don't search, or None - search.
      unpaired_msa, paired_msa, template_hits = _get_protein_msa_and_templates(
          sequence=chain.sequence,
          run_template_search=not has_templates,  # Skip template search if [].
          uniref90_msa_config=self._uniref90_msa_config,
          mgnify_msa_config=self._mgnify_msa_config,
          small_bfd_msa_config=self._small_bfd_msa_config,
          uniprot_msa_config=self._uniprot_msa_config,
          templates_config=self._templates_config,
          pdb_database_path=self._pdb_database_path,
      )
      unpaired_msa = unpaired_msa.to_a3m()
      paired_msa = paired_msa.to_a3m()
      templates = [
          folding_input.Template(
              mmcif=struc.to_mmcif(),
              query_to_template_map=hit.query_to_hit_mapping,
          )
          for hit, struc in template_hits.get_hits_with_structures()
      ]
    elif has_unpaired_msa and has_paired_msa and not has_templates:
      # Has MSA, but doesn't have templates. Search for templates only.
      empty_msa = msa.Msa.from_empty(
          query_sequence=chain.sequence,
          chain_poly_type=mmcif_names.PROTEIN_CHAIN,
      ).to_a3m()
      unpaired_msa = chain.unpaired_msa or empty_msa
      paired_msa = chain.paired_msa or empty_msa
      template_hits = _get_protein_templates(
          sequence=chain.sequence,
          input_msa_a3m=unpaired_msa,
          run_template_search=True,
          templates_config=self._templates_config,
          pdb_database_path=self._pdb_database_path,
      )
      templates = [
          folding_input.Template(
              mmcif=struc.to_mmcif(),
              query_to_template_map=hit.query_to_hit_mapping,
          )
          for hit, struc in template_hits.get_hits_with_structures()
      ]
    else:
      # Has MSA and templates, don't search for anything.
      if not has_unpaired_msa or not has_paired_msa or not has_templates:
        raise ValueError(
            f'Protein chain {chain.id} has unpaired MSA, paired MSA, or'
            ' templates set only partially. If you want to run the pipeline'
            ' with custom MSA/templates, you need to set all of them. You can'
            ' set MSA to empty string and templates to empty list to signify'
            ' that they should not be used and searched for.'
        )
      logging.info(
          'Skipping MSA and template search for protein chain %s because it '
          'already has MSAs and templates.',
          chain.id,
      )
      if not chain.unpaired_msa:
        logging.info('Using empty unpaired MSA for protein chain %s', chain.id)
      if not chain.paired_msa:
        logging.info('Using empty paired MSA for protein chain %s', chain.id)
      if not chain.templates:
        logging.info('Using no templates for protein chain %s', chain.id)
      empty_msa = msa.Msa.from_empty(
          query_sequence=chain.sequence,
          chain_poly_type=mmcif_names.PROTEIN_CHAIN,
      ).to_a3m()
      unpaired_msa = chain.unpaired_msa or empty_msa
      paired_msa = chain.paired_msa or empty_msa
      templates = chain.templates

    return folding_input.ProteinChain(
        id=chain.id,
        sequence=chain.sequence,
        ptms=chain.ptms,
        unpaired_msa=unpaired_msa,
        paired_msa=paired_msa,
        templates=templates,
    )

  def process_rna_chain(
      self, chain: folding_input.RnaChain
  ) -> folding_input.RnaChain:
    """Processes a single RNA chain."""
    if chain.unpaired_msa is not None:
      # Don't run MSA tools if the chain already has an MSA.
      logging.info(
          'Skipping MSA search for RNA chain %s because it already has MSA.',
          chain.id,
      )
      if not chain.unpaired_msa:
        logging.info('Using empty unpaired MSA for RNA chain %s', chain.id)
      empty_msa = msa.Msa.from_empty(
          query_sequence=chain.sequence, chain_poly_type=mmcif_names.RNA_CHAIN
      ).to_a3m()
      unpaired_msa = chain.unpaired_msa or empty_msa
    else:
      unpaired_msa = _get_rna_msa(
          sequence=chain.sequence,
          nt_rna_msa_config=self._nt_rna_msa_config,
          rfam_msa_config=self._rfam_msa_config,
          rnacentral_msa_config=self._rnacentral_msa_config,
      ).to_a3m()
    return folding_input.RnaChain(
        id=chain.id,
        sequence=chain.sequence,
        modifications=chain.modifications,
        unpaired_msa=unpaired_msa,
    )

  def process(self, fold_input: folding_input.Input) -> folding_input.Input:
    """Runs MSA and template tools and returns a new Input with the results."""
    processed_chains = []
    for chain in fold_input.chains:
      print(f'Running data pipeline for chain {chain.id}...')
      process_chain_start_time = time.time()
      match chain:
        case folding_input.ProteinChain():
          processed_chains.append(self.process_protein_chain(chain))
        case folding_input.RnaChain():
          processed_chains.append(self.process_rna_chain(chain))
        case _:
          processed_chains.append(chain)
      print(
          f'Running data pipeline for chain {chain.id} took'
          f' {time.time() - process_chain_start_time:.2f} seconds',
      )

    return dataclasses.replace(fold_input, chains=processed_chains)
