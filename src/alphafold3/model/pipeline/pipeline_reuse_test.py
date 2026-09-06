# SPDX-License-Identifier: Apache-2.0
"""Regression tests for request-local seed-invariant feature reuse."""
import dataclasses
import json
import pathlib
from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.data import featurisation
from alphafold3.data.tools import rdkit_utils
from alphafold3.model import msa_pairing
from alphafold3.model import features
from alphafold3.model.pipeline import pipeline
import numpy as np


class ReuseTest(parameterized.TestCase):

  def assertExact(self, actual, expected):
    self.assertIs(type(actual), type(expected))
    if isinstance(actual, np.ndarray):
      self.assertEqual(actual.dtype, expected.dtype)
      self.assertEqual(actual.shape, expected.shape)
      self.assertEqual(actual.strides, expected.strides)
      self.assertEqual(actual.flags.writeable, expected.flags.writeable)
      if actual.dtype.hasobject:
        for a, e in zip(actual.flat, expected.flat, strict=True):
          self.assertExact(a, e)
      else:
        np.testing.assert_array_equal(actual, expected)
    elif dataclasses.is_dataclass(actual):
      for f in dataclasses.fields(actual):
        self.assertExact(getattr(actual, f.name), getattr(expected, f.name))
    elif isinstance(actual, dict):
      self.assertEqual(actual.keys(), expected.keys())
      for k in actual:
        self.assertExact(actual[k], expected[k])
    elif isinstance(actual, (list, tuple)):
      self.assertLen(actual, len(expected))
      for a, e in zip(actual, expected, strict=True):
        self.assertExact(a, e)
    elif hasattr(actual, 'get_table'):
      for name in actual.tables:
        self.assertExact(actual.get_table(name), expected.get_table(name))
      for name in (
          'name', 'release_date', 'resolution', 'structure_method',
          'bioassembly_data', 'chemical_components_data',
      ):
        self.assertExact(getattr(actual, name), getattr(expected, name))
    else:
      self.assertEqual(actual, expected)

  def make_input(self, ligand, seeds):
    seq = 'ACDEFGHIKLMNPQRSTVWY'
    chains = [{'protein': {
        'id': 'A', 'sequence': seq,
        'unpairedMsa': f'>query\n{seq}\n>hit\nCCDEFGHIKLMNPQRSTVWY\n',
        'pairedMsa': '', 'templates': [],
    }}]
    if ligand:
      chains.extend([
          {'protein': {'id': 'B', 'sequence': 'ACDE',
                       'unpairedMsa': '', 'pairedMsa': '', 'templates': []}},
          {'ligand': {'id': 'C', 'smiles': 'CCO'}},
      ])
    return folding_input.Input.from_json(json.dumps({
        'name': 'reuse', 'dialect': 'alphafold3', 'version': 1,
        'sequences': chains, 'modelSeeds': seeds,
    }))

  @parameterized.product(ligand=(False, True), frames=(False, True), n=(1, 4))
  def test_exact_reuse_and_work_counts(self, ligand, frames, n):
    inp = self.make_input(ligand, [0, 12312837, 2**32 - 1, 1][:n])
    ccd = chemical_components.Ccd()
    p = pipeline.WholePdbPipeline(config=pipeline.WholePdbPipeline.Config(
        deterministic_frames=frames, buckets=[128],
    ))
    expected = [p.process_item(inp, np.random.RandomState(s), ccd, s)
                for s in inp.rng_seeds]
    with (
        mock.patch.object(features.MSA, 'compute_features',
                          wraps=features.MSA.compute_features) as msa,
        mock.patch.object(features.RefStructure, 'compute_features',
                          wraps=features.RefStructure.compute_features) as ref,
    ):
      actual = list(p._process_items(fold_input=inp, ccd=ccd))
    self.assertExact(actual, expected)
    self.assertEqual(msa.call_count, 1)
    self.assertEqual(ref.call_count, n * (2 if frames else 1))
    for i in range(1, n):
      for key in features.MSA.from_data_dict(actual[0]).as_data_dict():
        self.assertFalse(np.shares_memory(actual[0][key], actual[i][key]))
      self.assertFalse(np.shares_memory(
          actual[0]['frames_mask'], actual[i]['frames_mask']))
    # Another request using the same pipeline must start with fresh reuse state.
    other = self.make_input(not ligand, [1, 2])
    other_actual = list(p._process_items(fold_input=other, ccd=ccd))
    other_expected = [p.process_item(other, np.random.RandomState(s), ccd, s)
                      for s in other.rng_seeds]
    self.assertExact(other_actual, other_expected)




  def test_nan_validation_preserved(self):
    inp = self.make_input(False, [1, 2])
    ccd = chemical_components.Ccd()
    p = pipeline.WholePdbPipeline(config=pipeline.WholePdbPipeline.Config())
    compute = features.MSA.compute_features

    def invalid(*args, **kwargs):
      msa = compute(*args, **kwargs)
      return dataclasses.replace(msa, profile=np.full_like(msa.profile, np.nan))

    with mock.patch.object(features.MSA, 'compute_features', side_effect=invalid):
      with self.assertRaisesRegex(pipeline.NanDataError, 'NaN feature: profile'):
        list(p._process_items(fold_input=inp, ccd=ccd))

  @parameterized.product(paired=(False, True), frames=(False, True))
  def test_deep_and_paired_msa_reuse(self, paired, frames):
    data = json.loads(self.make_input(False, [1, 2]).to_json())
    seq = data['sequences'][0]['protein']['sequence']
    if paired:
      data['sequences'] = []
      for chain_id, sequence in (('A', seq), ('B', seq[::-1])):
        human = 'A' + sequence[1:]
        mouse = sequence[:-1] + 'A'
        msa = (f'>query\n{sequence}\n>sp|P12345|TEST_HUMAN\n{human}\n'
               f'>sp|P67890|TEST_MOUSE\n{mouse}\n')
        data['sequences'].append({'protein': {
            'id': chain_id, 'sequence': sequence, 'pairedMsa': msa,
            'unpairedMsa': f'>query\n{sequence}\n', 'templates': [],
        }})
    else:
      rows = [f'>query\n{seq}\n']
      for i in range(16383):
        prefix = ''.join(seq[(i // 20**j) % 20] for j in range(4))
        rows.append(f'>hit_{i}\n{prefix}{seq[4:]}\n')
      data['sequences'][0]['protein']['unpairedMsa'] = ''.join(rows)
    inp = folding_input.Input.from_json(json.dumps(data))
    ccd = chemical_components.Ccd()
    p = pipeline.WholePdbPipeline(config=pipeline.WholePdbPipeline.Config(
        deterministic_frames=frames, buckets=[128],
    ))
    expected = [p.process_item(inp, np.random.RandomState(s), ccd, s)
                for s in inp.rng_seeds]
    with (
        mock.patch.object(features.MSA, 'compute_features',
                          wraps=features.MSA.compute_features) as msa,
        mock.patch.object(msa_pairing, 'create_paired_features',
                          wraps=msa_pairing.create_paired_features) as pairing,
    ):
      actual = list(p._process_items(fold_input=inp, ccd=ccd))
    self.assertExact(actual, expected)
    self.assertEqual(msa.call_count, 1)
    if paired:
      self.assertEqual(pairing.call_count, 1)
      depth = int(actual[0]['num_alignments'])
      tokens = int(actual[0]['seq_length'])
      self.assertGreaterEqual(
          np.all(actual[0]['msa_mask'][:depth, :tokens], axis=1).sum(), 3)
    else:
      self.assertGreater(int(actual[0]['num_alignments']), 1000)

  def test_msa_mutation_during_iteration_and_single_seed_ownership(self):
    inp = self.make_input(False, [1, 2])
    ccd = chemical_components.Ccd()
    p = pipeline.WholePdbPipeline(config=pipeline.WholePdbPipeline.Config())
    iterator = p._process_items(fold_input=inp, ccd=ccd)
    first = next(iterator)
    for value in features.MSA.from_data_dict(first).as_data_dict().values():
      value.flat[0] = 1 if value.flat[0] == 0 else 0
    expected = p.process_item(inp, np.random.RandomState(2), ccd, 2)
    self.assertExact(next(iterator), expected)
    iterator.close()
    captured = []
    compute = features.MSA.compute_features

    def capture(*args, **kwargs):
      value = compute(*args, **kwargs)
      captured.append(value)
      return value

    with mock.patch.object(features.MSA, 'compute_features', side_effect=capture):
      single, = p._process_items(
          fold_input=dataclasses.replace(inp, rng_seeds=[1]), ccd=ccd)
    self.assertLen(captured, 1)
    for key, value in captured[0].as_data_dict().items():
      self.assertIs(single[key], value)

  def test_public_api_and_none_seed_rng(self):
    inp = self.make_input(False, [1, 2])
    ccd = chemical_components.Ccd()
    p = pipeline.WholePdbPipeline(config=pipeline.WholePdbPipeline.Config())
    actual = featurisation.featurise_input(inp, ccd, buckets=None)
    self.assertExact(actual, list(p._process_items(fold_input=inp, ccd=ccd)))
    actual_rng = np.random.RandomState(123)
    expected_rng = np.random.RandomState(123)
    seed = expected_rng.randint(2**31)
    expected = p.process_item(inp, expected_rng, ccd, random_seed=seed)
    actual = p.process_item(inp, actual_rng, ccd, random_seed=None)
    self.assertExact(actual, expected)
    self.assertExact(actual_rng.get_state(), expected_rng.get_state())

  @parameterized.parameters('unpairedMsa', 'pairedMsa', 'templates')
  def test_missing_feature_fields_still_fail(self, field):
    data = json.loads(self.make_input(False, [1, 2]).to_json())
    del data['sequences'][0]['protein'][field]
    inp = folding_input.Input.from_json(json.dumps(data))
    with self.assertRaisesRegex(ValueError, 'missing'):
      featurisation.featurise_input(inp, chemical_components.Ccd(), buckets=None)



if __name__ == '__main__':
  absltest.main()
