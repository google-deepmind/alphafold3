# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

from collections.abc import Sequence
import copy
import dataclasses
import datetime
import os

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import mmcif_names
from alphafold3.data import msa
from alphafold3.data import msa_config
from alphafold3.data import pipeline

import shutil

_JACKHMMER_BINARY_PATH = shutil.which('jackhmmer')
_NHMMER_BINARY_PATH = shutil.which('nhmmer')
_HMMALIGN_BINARY_PATH = shutil.which('hmmalign')
_HMMSEARCH_BINARY_PATH = shutil.which('hmmsearch')
_HMMBUILD_BINARY_PATH = shutil.which('hmmbuild')


def create_fasta_file(
    descriptions: Sequence[str],
    sequences: Sequence[str],
    path: str,
) -> None:
  with open(path, 'w') as f:
    for description, sequence in zip(descriptions, sequences, strict=True):
      f.write(f'>{description}\n{sequence}\n')


class PipelineTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()

    seqres_database_path = data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/'
        / 'pdb_seqres_2022_09_28__subsampled_1000.fasta'
    ).path()
    rfam_database_path = data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/'
        / 'rfam_14_4_clustered_rep_seq__subsampled_1000.fasta'
    ).path()
    pdb_database_path = data.Data(
        resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
    ).path()

    self.data_pipeline_config = pipeline.DataPipelineConfig(
        jackhmmer_binary_path=_JACKHMMER_BINARY_PATH,
        nhmmer_binary_path=_NHMMER_BINARY_PATH,
        hmmalign_binary_path=_HMMALIGN_BINARY_PATH,
        hmmsearch_binary_path=_HMMSEARCH_BINARY_PATH,
        hmmbuild_binary_path=_HMMBUILD_BINARY_PATH,
        small_bfd_database_path=seqres_database_path,
        mgnify_database_path=seqres_database_path,
        uniprot_cluster_annot_database_path=seqres_database_path,
        uniref90_database_path=seqres_database_path,
        ntrna_database_path=rfam_database_path,
        rfam_database_path=rfam_database_path,
        rna_central_database_path=rfam_database_path,
        pdb_database_path=pdb_database_path,
        seqres_database_path=seqres_database_path,
        max_template_date=datetime.date(2021, 9, 30),
    )

  @parameterized.product(
      unpaired_msa=['', '>foo\nTSTU\nT---\n'],
      paired_msa=['', '>foo\nTSTP\nT---\n'],
  )
  def test_process_protein_chain_custom_msa(self, unpaired_msa, paired_msa):
    data_pipeline = pipeline.DataPipeline(self.data_pipeline_config)
    protein_chain = folding_input.ProteinChain(
        id='A',
        sequence='TST',
        ptms=[],
        unpaired_msa=unpaired_msa,
        paired_msa=paired_msa,
        templates=[],
    )
    new_protein_chain = data_pipeline.process_protein_chain(protein_chain)
    expected_protein_chain = folding_input.ProteinChain(
        id='A',
        sequence='TST',
        ptms=[],
        unpaired_msa=unpaired_msa or '>Original query\nTST\n',
        paired_msa=paired_msa or '>Original query\nTST\n',
        templates=[],
    )
    self.assertEqual(new_protein_chain, expected_protein_chain)

  @parameterized.parameters(
      ('>foo\nTST', None, None),
      (None, '>foo\nTST', None),
      ('>foo\nTST', None, []),
      (None, '>foo\nTST', []),
      (
          None,
          None,
          [folding_input.Template(mmcif='tst', query_to_template_map={})],
      ),
  )
  def test_process_protein_chain_custom_msa_error(
      self,
      unpaired_msa: str | None,
      paired_msa: str | None,
      templates: Sequence[folding_input.Template] | None,
  ):
    data_pipeline = pipeline.DataPipeline(self.data_pipeline_config)

    protein_chain = folding_input.ProteinChain(
        id='A',
        sequence='TST',
        ptms=[],
        unpaired_msa=unpaired_msa,
        paired_msa=paired_msa,
        templates=templates,
    )

    with self.assertRaisesRegex(ValueError, 'set only partially'):
      data_pipeline.process_protein_chain(protein_chain)

  @parameterized.named_parameters(
      ('empty', ''), ('nonempty', '>query\nGGG\n>alignment1\nG-G\n')
  )
  def test_process_rna_chain_custom_msa(self, custom_msa):
    data_pipeline = pipeline.DataPipeline(self.data_pipeline_config)
    rna_chain = folding_input.RnaChain(
        id='A',
        sequence='GGG',
        modifications=[],
        unpaired_msa=custom_msa,
    )
    expected_rna_chain = folding_input.RnaChain(
        id='A',
        sequence='GGG',
        modifications=[],
        unpaired_msa=custom_msa or '>Original query\nGGG\n',
    )
    new_rna_chain = data_pipeline.process_rna_chain(rna_chain)
    self.assertEqual(new_rna_chain, expected_rna_chain)

  def test_protein_msa_deduplication(self):
    """Specifically checks logic for protein MSA deduplication."""

    # Define the databases.
    # We will use a single base sequence, and add various modifications to it.
    # The AAs added as modifications are chosen specifically so that Jackhmmer
    # classifies (at least some of them) as substitutions, rather than
    # insertions or deletions.
    base = 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'
    query1 = f'{base}MMM'
    query2 = f'{base}MMG'
    uniref90_descriptions = ['query', 'uniref90_0', 'uniref90_1', 'uniref90_2']
    uniref90_sequences = [query1, f'{base}AAA', f'{base}AAA', f'{base}AAA']
    small_bfd_descriptions = ['query', 'bfd_0', 'bfd_1', 'bfd_2']
    small_bfd_sequences = [query1, f'{base}AAA', f'{base}CCC', f'{base}FFF']
    mgnify_descriptions = ['query', 'mgnify_0', 'mgnify_1', 'mgnify_2']
    mgnify_sequences = [query1, f'{base}AAA', f'{base}III', f'{base}CCC']
    uniprot_descriptions = ['query', 'uniprot_0', 'uniprot_1', 'uniprot_2']
    uniprot_sequences = [query1, f'{base}AAA', f'{base}AAA', f'{base}AAA']

    # Write the databases to disk.
    database_dir = self.create_tempdir()
    uniref90_database_path = os.path.join(database_dir, 'uniref90.fasta')
    small_bfd_database_path = os.path.join(database_dir, 'small_bfd.fasta')
    mgnify_database_path = os.path.join(database_dir, 'mgnify.fasta')
    uniprot_database_path = os.path.join(database_dir, 'uniprot.fasta')
    create_fasta_file(
        descriptions=uniref90_descriptions,
        sequences=uniref90_sequences,
        path=uniref90_database_path,
    )
    create_fasta_file(
        descriptions=small_bfd_descriptions,
        sequences=small_bfd_sequences,
        path=small_bfd_database_path,
    )
    create_fasta_file(
        descriptions=mgnify_descriptions,
        sequences=mgnify_sequences,
        path=mgnify_database_path,
    )
    create_fasta_file(
        descriptions=uniprot_descriptions,
        sequences=uniprot_sequences,
        path=uniprot_database_path,
    )

    # Create a small version of pdb seqres so that Hmmsearch completes quickly.
    # For this test we don't care about what Hmmsearch does anyways.
    seqres_database_path = os.path.join(database_dir, 'pdb_seqres.fasta')
    with open(seqres_database_path, 'w') as f:
      f.write(
          """>101m_A mol:protein length:154  MYOGLOBIN
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
>102l_A mol:protein length:165  T4 LYSOZYME
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL"""
      )

    # Instantiate the pipeline, using the databases we just created.
    data_pipeline = pipeline.DataPipeline(
        dataclasses.replace(
            self.data_pipeline_config,
            uniref90_database_path=uniref90_database_path,
            uniprot_cluster_annot_database_path=uniprot_database_path,
            small_bfd_database_path=small_bfd_database_path,
            mgnify_database_path=mgnify_database_path,
            seqres_database_path=seqres_database_path,
        )
    )

    # Run the pipeline (in particular MSA and subsequent deduplication).
    data_pipeline_output = data_pipeline.process(
        folding_input.Input(
            name='test_folding_input',
            rng_seeds=[42],
            chains=[
                folding_input.ProteinChain(id='A', sequence=query1, ptms=[]),
                folding_input.ProteinChain(id='B', sequence=query1, ptms=[]),
                folding_input.ProteinChain(id='C', sequence=query1, ptms=[]),
                folding_input.ProteinChain(id='D', sequence=query2, ptms=[]),
            ],
        )
    )

    # Check the results.

    # Since each sequence description indicates which database it came
    # from, we can check the priority of deduplication across
    # databases is correct.
    query1_a3m = '\n'.join([
        '>query',
        query1,
        f'>uniref90_0/1-{len(base) + 1} [subseq from] uniref90_0',
        f'{base}A--',
        f'>bfd_1/1-{len(base) + 1} [subseq from] bfd_1',
        f'{base}C--',
        f'>bfd_2/1-{len(base) + 1} [subseq from] bfd_2',
        f'{base}F--',
        f'>mgnify_1/1-{len(base) + 2} [subseq from] mgnify_1',
        f'{base}II-',
    ])
    query2_a3m = '\n'.join([
        '>query',
        query2,
        f'>query/1-{len(base) + 2} [subseq from] query',
        f'{base}MM-',
        f'>uniref90_0/1-{len(base) + 1} [subseq from] uniref90_0',
        f'{base}A--',
        f'>bfd_1/1-{len(base) + 1} [subseq from] bfd_1',
        f'{base}C--',
        f'>bfd_2/1-{len(base) + 1} [subseq from] bfd_2',
        f'{base}F--',
        f'>mgnify_1/1-{len(base) + 1} [subseq from] mgnify_1',
        f'{base}I--',
    ])

    chains = data_pipeline_output.protein_chains
    assert chains[0].unpaired_msa is not None
    assert chains[1].unpaired_msa is not None
    assert chains[2].unpaired_msa is not None
    assert chains[3].unpaired_msa is not None
    self.assertEqual(chains[0].unpaired_msa.strip(), query1_a3m)
    self.assertEqual(chains[1].unpaired_msa.strip(), query1_a3m)
    self.assertEqual(chains[2].unpaired_msa.strip(), query1_a3m)
    self.assertEqual(chains[3].unpaired_msa.strip(), query2_a3m)

  def test_rna_msa_deduplication(self):
    """Specifically checks logic for RNA MSA deduplication."""

    # Define the databases.
    # We will use a single base sequence, and add various modifications to it.
    base = 'GUGCAUGCCUAUAGCAUAUCAUUAAUGCACCAGAUCCCAUUAUAACUCCUCAUGUAAGCGUGCUCGAGAUAGAUUAGUACUGGGAUGGUUGACUGCAAAGGAAGUCUUAGUGUUUUAC'
    query = f'{base}AUG'
    rfam_descs = ['query', 'rfam_0', 'rfam_1', 'rfam_2']
    rfam_seqs = [query, f'{base}ACG', f'{base}ACG', f'{base}ACG']
    rnacentral_descs = ['query', 'rnacentral_0', 'rnacentral_1', 'rnacentral_2']
    rnacentral_seqs = [query, f'{base}ACG', f'{base}UCG', f'{base}UCG']
    ntrna_descs = ['query', 'ntrna_0', 'ntrna_1', 'ntrna_2']
    ntrna_seqs = [query, f'{base}ACG', f'{base}UCG', f'{base}AUA']

    # Write the databases to disk.
    database_dir = self.create_tempdir()
    rfam_db_path = os.path.join(database_dir, 'rfam.fasta')
    rnacentral_db_path = os.path.join(database_dir, 'rnacentral.fasta')
    ntrna_db_path = os.path.join(database_dir, 'ntrna.fasta')
    create_fasta_file(
        descriptions=rfam_descs, sequences=rfam_seqs, path=rfam_db_path
    )
    create_fasta_file(
        descriptions=rnacentral_descs,
        sequences=rnacentral_seqs,
        path=rnacentral_db_path,
    )
    create_fasta_file(
        descriptions=ntrna_descs, sequences=ntrna_seqs, path=ntrna_db_path
    )

    # Instantiate the pipeline, using the databases we just created.
    data_pipeline = pipeline.DataPipeline(
        dataclasses.replace(
            self.data_pipeline_config,
            ntrna_database_path=ntrna_db_path,
            rfam_database_path=rfam_db_path,
            rna_central_database_path=rnacentral_db_path,
        )
    )

    # Run the pipeline (in particular MSA and subsequent deduplication).
    new_rna_chain = data_pipeline.process_rna_chain(
        folding_input.RnaChain(id='test', sequence=query, modifications=[])
    )

    # Check the results.
    # Since each sequence description indicates which database it came from,
    # we can check that the priority of deduplication across databases is
    # correct.
    expected_a3m = '\n'.join([
        '>query',
        query,
        f'>rfam_0/1-{len(base) + 2} [subseq from] rfam_0',
        f'{base}A--c',
        f'>rnacentral_1/1-{len(base)} [subseq from] rnacentral_1',
        f'{base}---',
        f'>ntrna_2/1-{len(base) + 2} [subseq from] ntrna_2',
        f'{base}AU-',
    ])
    assert new_rna_chain.unpaired_msa is not None
    self.assertEqual(expected_a3m, new_rna_chain.unpaired_msa.strip())


if __name__ == '__main__':
  absltest.main()
