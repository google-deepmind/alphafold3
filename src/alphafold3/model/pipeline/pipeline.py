# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""The main featurizer."""

import bisect
from collections.abc import Sequence
import datetime
import itertools

from absl import logging
from alphafold3.common import base_config
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.model import feat_batch
from alphafold3.model import features
from alphafold3.model.pipeline import inter_chain_bonds
from alphafold3.model.pipeline import structure_cleaning
from alphafold3.structure import chemical_components as struc_chem_comps
import numpy as np


_DETERMINISTIC_FRAMES_RANDOM_SEED = 12312837


def calculate_bucket_size(
    num_tokens: int, buckets: Sequence[int] | None
) -> int:
  """Calculates the bucket size to pad the data to."""
  if buckets is None:
    return num_tokens

  if not buckets:
    raise ValueError('Buckets must be non-empty.')

  if not all(prev < curr for prev, curr in itertools.pairwise(buckets)):
    raise ValueError(
        f'Buckets must be in strictly increasing order. Got {buckets=}.'
    )

  bucket_idx = bisect.bisect_left(buckets, num_tokens)

  if bucket_idx == len(buckets):
    logging.warning(
        'Creating a new bucket of size %d since the input has more tokens than'
        ' the largest bucket size %d. This may trigger a re-compilation of the'
        ' model. Consider additional large bucket sizes to avoid excessive'
        ' re-compilation.',
        num_tokens,
        buckets[-1],
    )
    return num_tokens

  return buckets[bucket_idx]


class NanDataError(Exception):
  """Raised if the data pipeline produces data containing nans."""


class TotalNumResOutOfRangeError(Exception):
  """Raised if total number of residues for all chains outside allowed range."""


class MmcifNumChainsError(Exception):
  """Raised if the mmcif file contains too many / too few chains."""


class WholePdbPipeline:
  """Processes an entire mmcif entity and merges the content."""

  class Config(base_config.BaseConfig):
    """Configuration object for `WholePdbPipeline`.

    Properties:
      max_atoms_per_token: number of atom slots in one token (was called
        num_dense, and semi-hardcoded to 24 before)
      pad_num_chains: Size to pad NUM_CHAINS feature dimensions to, only for
        protein chains.
      buckets: Bucket sizes to pad the data to, to avoid excessive
        re-compilation of the model. If None, calculate the appropriate bucket
        size from the number of tokens. If not None, must be a sequence of at
        least one integer, in strictly increasing order. Will raise an error if
        the number of tokens is more than the largest bucket size.
      max_total_residues: Any mmCIF with more total residues will be rejected.
        If none, then no limit is applied.
      min_total_residues: Any mmCIF with less total residues will be rejected.
      msa_crop_size: Maximum size of MSA to take across all chains.
      max_template_date: Optional max template date to prevent data leakage in
        validation.
      ref_max_modified_date: Optional maximum date that controls whether to
        allow use of model coordinates for a chemical component from the CCD if
        RDKit conformer generation fails and the component does not have ideal
        coordinates set. Only for components that have been released before this
        date the model coordinates can be used as a fallback.
      max_templates: The maximum number of templates to send through the network
        set to 0 to switch off templates.
      filter_clashes: If true then will remove clashing chains.
      filter_crystal_aids: If true ligands in the cryal aid list are removed.
      max_paired_sequence_per_species: The maximum number of sequences per
        species that will be used for MSA pairing.
      drop_ligand_leaving_atoms: Flag for handling leaving atoms for ligands.
      average_num_atoms_per_token: Target average number of atoms per token to
        compute the padding size for flat atoms.
      atom_cross_att_queries_subset_size: queries subset size in atom cross
        attention
      atom_cross_att_keys_subset_size: keys subset size in atom cross attention
      flatten_non_standard_residues: Whether to expand non-standard polymer
        residues into flat-atom format.
      remove_nonsymmetric_bonds: Whether to remove nonsymmetric bonds from
        symmetric polymer chains.
      deterministic_frames: Whether to use fixed-seed reference positions to
        construct deterministic frames.
    """

    max_atoms_per_token: int = 24
    pad_num_chains: int = 1000
    buckets: list[int] | None = None
    max_total_residues: int | None = None
    min_total_residues: int | None = None
    msa_crop_size: int = 16384
    max_template_date: datetime.date | None = None
    ref_max_modified_date: datetime.date | None = None
    max_templates: int = 4
    filter_clashes: bool = False
    filter_crystal_aids: bool = False
    max_paired_sequence_per_species: int = 600
    drop_ligand_leaving_atoms: bool = True
    average_num_atoms_per_token: int = 24
    atom_cross_att_queries_subset_size: int = 32
    atom_cross_att_keys_subset_size: int = 128
    flatten_non_standard_residues: bool = True
    remove_nonsymmetric_bonds: bool = False
    deterministic_frames: bool = True
    conformer_max_iterations: int | None = None

  def __init__(self, *, config: Config):
    """Initializes WholePdb data pipeline.

    Args:
      config: Pipeline configuration.
    """
    self._config = config

  def process_item(
      self,
      fold_input: folding_input.Input,
      random_state: np.random.RandomState,
      ccd: chemical_components.Ccd,
      random_seed: int | None = None,
  ) -> features.BatchDict:
    """Takes requests from in_queue, adds (key, serialized ex) to out_queue."""
    if random_seed is None:
      random_seed = random_state.randint(2**31)

    random_state = np.random.RandomState(seed=random_seed)

    logging_name = f'{fold_input.name}, random_seed={random_seed}'
    logging.info('processing %s', logging_name)
    struct = fold_input.to_structure(ccd=ccd)

    # Clean structure.
    cleaned_struc, cleaning_metadata = structure_cleaning.clean_structure(
        struct,
        ccd=ccd,
        drop_non_standard_atoms=True,
        drop_missing_sequence=True,
        filter_clashes=self._config.filter_clashes,
        filter_crystal_aids=self._config.filter_crystal_aids,
        filter_waters=True,
        filter_hydrogens=True,
        filter_leaving_atoms=self._config.drop_ligand_leaving_atoms,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=True,
        remove_polymer_polymer_bonds=True,
        remove_bad_bonds=True,
        remove_nonsymmetric_bonds=self._config.remove_nonsymmetric_bonds,
    )

    num_clashing_chains_removed = cleaning_metadata[
        'num_clashing_chains_removed'
    ]

    if num_clashing_chains_removed:
      logging.info(
          'Removed %d clashing chains from %s',
          num_clashing_chains_removed,
          logging_name,
      )

    # No chains after fixes
    if cleaned_struc.num_chains == 0:
      raise MmcifNumChainsError(f'{logging_name}: No chains in structure!')

    polymer_ligand_bonds, ligand_ligand_bonds = (
        inter_chain_bonds.get_polymer_ligand_and_ligand_ligand_bonds(
            cleaned_struc,
            only_glycan_ligands=False,
            allow_multiple_bonds_per_atom=True,
        )
    )

    # If empty replace with None as this causes errors downstream.
    if ligand_ligand_bonds and not ligand_ligand_bonds.atom_name.size:
      ligand_ligand_bonds = None
    if polymer_ligand_bonds and not polymer_ligand_bonds.atom_name.size:
      polymer_ligand_bonds = None

    # Create the flat output AtomLayout
    empty_output_struc, flat_output_layout = (
        structure_cleaning.create_empty_output_struc_and_layout(
            struc=cleaned_struc,
            ccd=ccd,
            polymer_ligand_bonds=polymer_ligand_bonds,
            ligand_ligand_bonds=ligand_ligand_bonds,
            drop_ligand_leaving_atoms=self._config.drop_ligand_leaving_atoms,
        )
    )

    # Select the tokens for Evoformer.
    # Each token (e.g. a residue) is encoded as one representative atom. This
    # is flexible enough to allow the 1-token-per-atom ligand representation
    # in the future.
    all_tokens, all_token_atoms_layout, standard_token_idxs = (
        features.tokenizer(
            flat_output_layout,
            ccd=ccd,
            max_atoms_per_token=self._config.max_atoms_per_token,
            flatten_non_standard_residues=self._config.flatten_non_standard_residues,
            logging_name=logging_name,
        )
    )
    total_tokens = len(all_tokens.atom_name)
    if (
        self._config.max_total_residues
        and total_tokens > self._config.max_total_residues
    ):
      raise TotalNumResOutOfRangeError(
          'Total Number of Residues > max_total_residues: '
          f'({total_tokens} > {self._config.max_total_residues})'
      )

    if (
        self._config.min_total_residues
        and total_tokens < self._config.min_total_residues
    ):
      raise TotalNumResOutOfRangeError(
          'Total Number of Residues < min_total_residues: '
          f'({total_tokens} < {self._config.min_total_residues})'
      )

    logging.info(
        'Calculating bucket size for input with %d tokens.', total_tokens
    )
    padded_token_length = calculate_bucket_size(
        total_tokens, self._config.buckets
    )
    logging.info(
        'Got bucket size %d for input with %d tokens, resulting in %d padded'
        ' tokens.',
        padded_token_length,
        total_tokens,
        padded_token_length - total_tokens,
    )

    # Padding shapes for all features.
    num_atoms = padded_token_length * self._config.average_num_atoms_per_token
    # Round up to next multiple of subset size.
    num_atoms = int(
        np.ceil(num_atoms / self._config.atom_cross_att_queries_subset_size)
        * self._config.atom_cross_att_queries_subset_size
    )
    padding_shapes = features.PaddingShapes(
        num_tokens=padded_token_length,
        msa_size=self._config.msa_crop_size,
        num_chains=self._config.pad_num_chains,
        num_templates=self._config.max_templates,
        num_atoms=num_atoms,
    )

    # Create the atom layouts for flat atom cross attention
    batch_atom_cross_att = features.AtomCrossAtt.compute_features(
        all_token_atoms_layout=all_token_atoms_layout,
        queries_subset_size=self._config.atom_cross_att_queries_subset_size,
        keys_subset_size=self._config.atom_cross_att_keys_subset_size,
        padding_shapes=padding_shapes,
    )

    # Extract per-token features
    batch_token_features = features.TokenFeatures.compute_features(
        all_tokens=all_tokens,
        padding_shapes=padding_shapes,
    )

    # Create reference structure features
    chemical_components_data = struc_chem_comps.populate_missing_ccd_data(
        ccd=ccd,
        chemical_components_data=cleaned_struc.chemical_components_data,
        populate_pdbx_smiles=True,
    )

    # Add smiles info to empty_output_struc.
    empty_output_struc = empty_output_struc.copy_and_update_globals(
        chemical_components_data=chemical_components_data
    )
    # Create layouts and store structures for model output conversion.
    batch_convert_model_output = features.ConvertModelOutput.compute_features(
        all_token_atoms_layout=all_token_atoms_layout,
        padding_shapes=padding_shapes,
        cleaned_struc=cleaned_struc,
        flat_output_layout=flat_output_layout,
        empty_output_struc=empty_output_struc,
        polymer_ligand_bonds=polymer_ligand_bonds,
        ligand_ligand_bonds=ligand_ligand_bonds,
    )

    # Create the PredictedStructureInfo
    batch_predicted_structure_info = (
        features.PredictedStructureInfo.compute_features(
            all_tokens=all_tokens,
            all_token_atoms_layout=all_token_atoms_layout,
            padding_shapes=padding_shapes,
        )
    )

    # Create MSA features
    batch_msa = features.MSA.compute_features(
        all_tokens=all_tokens,
        standard_token_idxs=standard_token_idxs,
        padding_shapes=padding_shapes,
        fold_input=fold_input,
        logging_name=logging_name,
        max_paired_sequence_per_species=self._config.max_paired_sequence_per_species,
    )

    # Create template features
    batch_templates = features.Templates.compute_features(
        all_tokens=all_tokens,
        standard_token_idxs=standard_token_idxs,
        padding_shapes=padding_shapes,
        fold_input=fold_input,
        max_templates=self._config.max_templates,
        logging_name=logging_name,
    )

    ref_max_modified_date = self._config.ref_max_modified_date
    conformer_max_iterations = self._config.conformer_max_iterations
    batch_ref_structure, ligand_ligand_bonds = (
        features.RefStructure.compute_features(
            all_token_atoms_layout=all_token_atoms_layout,
            ccd=ccd,
            padding_shapes=padding_shapes,
            chemical_components_data=chemical_components_data,
            random_state=random_state,
            ref_max_modified_date=ref_max_modified_date,
            conformer_max_iterations=conformer_max_iterations,
            ligand_ligand_bonds=ligand_ligand_bonds,
        )
    )
    deterministic_ref_structure = None
    if self._config.deterministic_frames:
      deterministic_ref_structure, _ = features.RefStructure.compute_features(
          all_token_atoms_layout=all_token_atoms_layout,
          ccd=ccd,
          padding_shapes=padding_shapes,
          chemical_components_data=chemical_components_data,
          random_state=(
              np.random.RandomState(_DETERMINISTIC_FRAMES_RANDOM_SEED)
          ),
          ref_max_modified_date=ref_max_modified_date,
          conformer_max_iterations=None,
          ligand_ligand_bonds=ligand_ligand_bonds,
      )

    # Create ligand-polymer bond features.
    polymer_ligand_bond_info = features.PolymerLigandBondInfo.compute_features(
        all_tokens=all_tokens,
        all_token_atoms_layout=all_token_atoms_layout,
        bond_layout=polymer_ligand_bonds,
        padding_shapes=padding_shapes,
    )
    # Create ligand-ligand bond features.
    ligand_ligand_bond_info = features.LigandLigandBondInfo.compute_features(
        all_tokens,
        ligand_ligand_bonds,
        padding_shapes,
    )

    # Create the Pseudo-beta layout for distogram head and distance error head.
    batch_pseudo_beta_info = features.PseudoBetaInfo.compute_features(
        all_token_atoms_layout=all_token_atoms_layout,
        ccd=ccd,
        padding_shapes=padding_shapes,
        logging_name=logging_name,
    )

    # Frame construction.
    batch_frames = features.Frames.compute_features(
        all_tokens=all_tokens,
        all_token_atoms_layout=all_token_atoms_layout,
        ref_structure=(
            deterministic_ref_structure
            if self._config.deterministic_frames
            else batch_ref_structure
        ),
        padding_shapes=padding_shapes,
    )

    # Assemble the Batch object.
    batch = feat_batch.Batch(
        msa=batch_msa,
        templates=batch_templates,
        token_features=batch_token_features,
        ref_structure=batch_ref_structure,
        predicted_structure_info=batch_predicted_structure_info,
        polymer_ligand_bond_info=polymer_ligand_bond_info,
        ligand_ligand_bond_info=ligand_ligand_bond_info,
        pseudo_beta_info=batch_pseudo_beta_info,
        atom_cross_att=batch_atom_cross_att,
        convert_model_output=batch_convert_model_output,
        frames=batch_frames,
    )

    np_example = batch.as_data_dict()
    if 'num_iter_recycling' in np_example:
      del np_example['num_iter_recycling']  # that does not belong here

    for name, value in np_example.items():
      if (
          value.dtype.kind not in {'U', 'S'}
          and value.dtype.name != 'object'
          and np.isnan(np.sum(value))
      ):
        raise NanDataError(
            'The output of the data pipeline contained nans. '
            f'nan feature: {name}, fold input name: {fold_input.name}, '
            f'random_seed {random_seed}'
        )

    return np_example
