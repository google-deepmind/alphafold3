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

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.model import confidences
import jax.numpy as jnp
import numpy as np


_TEST_DATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)


class ConfidencesTest(parameterized.TestCase):

  @parameterized.parameters(
      ('2pv7', 0.0),
      ('2fid', 0.13533834586466165),
      ('6gnx', 0.24427480916030533),
      ('1wib', 0.41304347826086957),
      ('4bxt', 0.6666666666666666),
      ('1al4', 1.0),
  )
  def test_fraction_disordered(self, cif_name, exp_fraction_disordered):
    cif_str = _TEST_DATA.load(f'{cif_name}.cif')
    struc = structure.from_mmcif(cif_str)
    fraction_disordered = confidences.fraction_disordered(struc=struc)
    self.assertAlmostEqual(fraction_disordered, exp_fraction_disordered)

  def test_fraction_disordered_with_multi_letter_chain_ids(self):
    cif_str = _TEST_DATA.load('3b9f.cif')
    struc = structure.from_mmcif(cif_str)
    struc = struc.rename_chain_ids(new_id_by_old_id={'A': 'AA', 'B': 'BB'})
    fraction_disordered = confidences.fraction_disordered(struc=struc)
    self.assertAlmostEqual(fraction_disordered, 0.05666156202143951)

  def test_ranking_metric(self):
    full_pde = jnp.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
    contact_probs = jnp.array([[0.5, 1.0], [0.5, 0.0]])
    rank_metric = confidences.rank_metric(
        full_pde=full_pde, contact_probs=contact_probs
    )
    np.testing.assert_allclose(rank_metric, np.array([-2.0, -6.0]), atol=1e-3)

  def test_chain_pair_pde(self):
    full_pde = np.array(
        [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[2, 3, 4], [5, 6, 7], [8, 9, 10]]]
    )
    asym_ids = np.array([1, 1, 2])
    pde_mean, pde_min = confidences.chain_pair_pde(
        num_tokens=3, asym_ids=asym_ids, full_pde=full_pde
    )
    np.testing.assert_allclose(
        pde_mean, np.array([[[3, 4.5], [7.5, 9]], [[4, 5.5], [8.5, 10.0]]])
    )
    np.testing.assert_allclose(
        pde_min, np.array([[[1.0, 3.0], [7.0, 9.0]], [[2.0, 4.0], [8.0, 10.0]]])
    )

  def test_pde_single(self):
    full_pde = np.array(
        [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[2, 3, 4], [5, 6, 7], [8, 9, 10]]]
    )
    contact_probs = np.array([[0, 1, 1], [1, 1, 0], [0, 1, 0]])
    asym_ids = np.array([1, 1, 2])
    (intra_chain_single_pde, cross_chain_single_pde, all_tokens_single_pde) = (
        confidences.pde_single(
            num_tokens=3,
            asym_ids=asym_ids,
            full_pde=full_pde,
            contact_probs=contact_probs,
        )
    )
    np.testing.assert_allclose(
        intra_chain_single_pde, np.array([[3.666666666, 0], [4.666666666, 0]])
    )
    np.testing.assert_allclose(
        cross_chain_single_pde, np.array([[3, 8], [4, 9]])
    )
    np.testing.assert_allclose(
        all_tokens_single_pde, np.array([[2.5, 4.5, 8], [3.5, 5.5, 9]])
    )

  def test_chain_pair_pae_mask(self):
    full_pae = np.array(
        [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[2, 3, 4], [5, 6, 7], [8, 9, 10]]]
    )
    mask = np.array([[0, 0, 0], [1, 1, 0], [1, 0, 1]]).astype(bool)
    asym_ids = np.array([1, 1, 2])
    pe_mean, pe_min, _ = confidences.chain_pair_pae(
        num_tokens=3, asym_ids=asym_ids, full_pae=full_pae, mask=mask
    )
    np.testing.assert_allclose(
        pe_mean,
        np.array([[[4.5, np.nan], [7.0, 9.0]], [[5.5, np.nan], [8.0, 10.0]]]),
    )
    np.testing.assert_allclose(
        pe_min, np.array([[[4, np.nan], [7, 9]], [[5, np.nan], [8, 10]]])
    )

  def test_chain_pair_pae_contact_probs(self):
    full_pae = np.array([
        [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ],
        [
            [2, 3, 4],
            [5, 6, 7],
            [8, 9, 10],
        ],
    ])
    mask = np.array([
        [0, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
    ]).astype(bool)
    contact_probs = np.array([
        [10, 10, 10],
        [1, 2, 0],
        [1, 10, 2],
    ])
    asym_ids = np.array([1, 1, 2])
    pe_mean, pe_min, _ = confidences.chain_pair_pae(
        num_tokens=3,
        asym_ids=asym_ids,
        full_pae=full_pae,
        mask=mask,
        contact_probs=contact_probs,
    )
    np.testing.assert_allclose(
        pe_mean,
        np.array([
            [[4.666667, np.nan], [7.0, 9.0]],
            [[5.666667, np.nan], [8.0, 10.0]],
        ]),
    )
    np.testing.assert_allclose(
        pe_min, np.array([[[4, np.nan], [7, 9]], [[5, np.nan], [8, 10]]])
    )

  def test_weighted_nanmean(self):
    value = np.array([[1, 1], [2, 2], [np.nan, np.nan]], dtype=float)
    mask = np.array([[1, 1], [1 / 2, 1 / 2], [0.1, 0.1]], dtype=float)

    np.testing.assert_allclose(
        confidences.weighted_nanmean(value, mask, axis=0),
        [1.3333333, 1.3333333],
    )
    np.testing.assert_allclose(
        confidences.weighted_nanmean(value, mask, axis=1), [1.0, 2.0, np.nan]
    )

  @parameterized.parameters(
      ([np.nan, np.nan], [1, 1], 0, np.nan),
      ([[np.nan, np.nan]], [[0, 0]], 0, [np.nan, np.nan]),
      ([[np.nan, np.nan]], [[0, 0]], 1, [np.nan]),
      ([[0, np.nan]], [[np.nan, 0]], 0, [np.nan, np.nan]),
      ([[0, np.nan]], [[np.nan, 0]], 1, [np.nan]),
  )
  def test_weighted_nanmean_all_nan(self, value, mask, axis, expected):
    mean = confidences.weighted_nanmean(
        np.array(value, dtype=np.float64),
        np.array(mask, dtype=np.float64),
        axis=axis,
    )
    np.testing.assert_array_equal(
        mean, np.array(expected, dtype=np.float64), strict=True
    )

  def test_reduce_chain_pair(self):
    full_pae = np.array([
        [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ],
        [
            [2, 3, 4],
            [5, 6, 7],
            [8, 9, 10],
        ],
    ])
    mask = np.array([
        [0, 1, 1],
        [1, 1, 0],
        [0, 1, 0],
    ])
    asym_ids = np.array([1, 1, 2])
    chain_pair_pae, _, unique_asym_ids = confidences.chain_pair_pae(
        num_tokens=3,
        asym_ids=asym_ids,
        full_pae=full_pae,
        mask=mask,
        contact_probs=np.ones_like(mask),
    )
    num_chain_tokens = np.array(
        [sum(asym_ids == asym_id) for asym_id in unique_asym_ids]
    )
    # aggregate over columns
    (ichain_pae, xchain_pae) = confidences.reduce_chain_pair(
        chain_pair_met=chain_pair_pae,
        num_chain_tokens=num_chain_tokens,
        agg_over_col=True,
        agg_type='mean',
        weight_method='per_token',
    )
    np.testing.assert_allclose(
        ichain_pae, np.array([[3.666666666, np.nan], [4.666666666, np.nan]])
    )
    np.testing.assert_allclose(
        xchain_pae,
        np.array([
            [3, 8],
            [4, 9],
        ]),
    )
    # aggregate over rows
    (ichain_pae, xchain_pae) = confidences.reduce_chain_pair(
        chain_pair_met=chain_pair_pae,
        num_chain_tokens=num_chain_tokens,
        agg_over_col=False,
        agg_type='mean',
        weight_method='per_token',
    )
    np.testing.assert_allclose(
        ichain_pae, np.array([[3.666666666, np.nan], [4.666666666, np.nan]])
    )
    np.testing.assert_allclose(xchain_pae, np.array([[8, 3], [9, 4]]))

    # aggregate over rows, min
    full_pae = np.array([
        [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ],
        [
            [2, 3, 4],
            [5, 6, 7],
            [8, 9, 10],
        ],
    ])
    mask = np.array([
        [0, 1, 1],
        [1, 1, 0],
        [0, 0, 0],
    ])
    asym_ids = np.array([1, 2, 3])
    chain_pair_pae, _, unique_asym_ids = confidences.chain_pair_pae(
        num_tokens=3,
        asym_ids=asym_ids,
        full_pae=full_pae,
        mask=mask,
        contact_probs=np.ones_like(mask),
    )
    num_chain_tokens = np.array(
        [sum(asym_ids == asym_id) for asym_id in unique_asym_ids]
    )
    (ichain_pae, xchain_pae) = confidences.reduce_chain_pair(
        chain_pair_met=chain_pair_pae,
        num_chain_tokens=num_chain_tokens,
        agg_over_col=False,
        agg_type='min',
        weight_method='per_token',
    )
    np.testing.assert_allclose(
        ichain_pae, np.array([[np.nan, 5, np.nan], [np.nan, 6, np.nan]])
    )
    np.testing.assert_allclose(xchain_pae, np.array([[4, 2, 3], [5, 3, 4]]))
    # aggregate over cols, min
    (ichain_pae, xchain_pae) = confidences.reduce_chain_pair(
        chain_pair_met=chain_pair_pae,
        num_chain_tokens=num_chain_tokens,
        agg_over_col=True,
        agg_type='min',
        weight_method='per_token',
    )
    np.testing.assert_allclose(
        ichain_pae, np.array([[np.nan, 5, np.nan], [np.nan, 6, np.nan]])
    )
    np.testing.assert_allclose(
        xchain_pae, np.array([[2, 4, np.nan], [3, 5, np.nan]])
    )

    full_pae = np.array([
        [
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            [2, 2, 2, 2],
            [3, 3, 3, 3],
        ],
    ])
    mask = np.ones_like(full_pae, shape=full_pae.shape[1:])
    asym_ids = np.array([1, 1, 2, 3])
    chain_pair_pae, _, unique_asym_ids = confidences.chain_pair_pae(
        num_tokens=full_pae.shape[-1],
        asym_ids=asym_ids,
        full_pae=full_pae,
        mask=mask,
        contact_probs=np.ones_like(mask),
    )
    num_chain_tokens = np.array(
        [sum(asym_ids == asym_id) for asym_id in unique_asym_ids]
    )

    def fn(agg_over_col: bool, weight_method: str):
      (
          _,
          ret,
      ) = confidences.reduce_chain_pair(
          chain_pair_met=chain_pair_pae,
          num_chain_tokens=num_chain_tokens,
          agg_over_col=agg_over_col,
          agg_type='mean',
          weight_method=weight_method,
      )
      return ret

    np.testing.assert_allclose(
        fn(agg_over_col=False, weight_method='per_chain'),
        np.array([[2.5, 2, 1.5]]),
    )
    np.testing.assert_allclose(
        fn(agg_over_col=True, weight_method='per_chain'),
        np.array([[1, 2, 3]]),
    )
    np.testing.assert_allclose(
        fn(agg_over_col=False, weight_method='per_token'),
        np.array([[2.5, 1.6666666, 1.3333333]]),
    )
    np.testing.assert_allclose(
        fn(agg_over_col=True, weight_method='per_token'),
        np.array([[1, 2, 3]]),
    )


if __name__ == '__main__':
  absltest.main()
