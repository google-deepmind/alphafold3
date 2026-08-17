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

"""Alignment based metrics."""

import numpy as np


def transform_ls(
    x: np.ndarray,
    b: np.ndarray,
    *,
    allow_reflection: bool = False,
) -> np.ndarray:
  """Find the least squares best fit rotation between two sets of N points.

  Solve Ax = b for A. Where A is the transform rotating x^T into b^T.

  Args:
    x: NxD numpy array of coordinates. Usually dimension D is 3.
    b: NxD numpy array of coordinates. Usually dimension D is 3.
    allow_reflection: Whether the returned transformation can reflect as well as
      rotate.

  Returns:
    Matrix A transforming x into b, i.e. s.t. Ax^T = b^T.
  """
  assert x.shape[1] >= b.shape[1]
  assert b.shape[0] == x.shape[0], '%d, %d' % (b.shape[0], x.shape[0])
  # First postmultiply by x.;
  # Axx^t = b x^t
  bxt = np.dot(b.transpose(), x) / b.shape[0]

  u, _, v = np.linalg.svd(bxt)

  r = np.dot(u, v)
  if not allow_reflection:
    flip = np.ones((v.shape[1], 1))
    flip[v.shape[1] - 1, 0] = np.sign(np.linalg.det(r))
    r = np.dot(u, v * flip)

  return r


def align(
    *,
    x: np.ndarray,
    y: np.ndarray,
    x_indices: np.ndarray,
    y_indices: np.ndarray,
) -> np.ndarray:
  """Align x to y considering only included_idxs.

  Args:
    x: NxD np array of coordinates.
    y: NxD np array of coordinates.
    x_indices: An np array of indices for `x` that will be used in the
      alignment. Must be of the same length as `y_included_idxs`.
    y_indices: An np array of indices for `y` that will be used in the
      alignment. Must be of the same length as `x_included_idxs`.

  Returns:
    NxD np array of points obtained by applying a rigid transformation to x.
    These points are aligned to y and the alignment is the optimal alignment
    over the points in included_idxs.

  Raises:
    ValueError: If the number of included indices is not the same for both
    input arrays.
  """
  if len(x_indices) != len(y_indices):
    raise ValueError(
        'Number of included indices must be the same for both input arrays,'
        f' but got for x: {len(x_indices)}, and for y: {len(y_indices)}.'
    )

  x_mean = np.mean(x[x_indices, :], axis=0)
  y_mean = np.mean(y[y_indices, :], axis=0)

  centered_x = x - x_mean
  centered_y = y - y_mean

  t = transform_ls(centered_x[x_indices, :], centered_y[y_indices, :])
  transformed_x = np.dot(centered_x, t.transpose()) + y_mean

  return transformed_x


def deviations_from_coords(
    decoy_coords: np.ndarray,
    gt_coords: np.ndarray,
    align_idxs: np.ndarray | None = None,
    include_idxs: np.ndarray | None = None,
) -> np.ndarray:
  """Returns the raw per-atom deviations used in RMSD computation."""
  if decoy_coords.shape != gt_coords.shape:
    raise ValueError(
        'decoy_coords.shape and gt_coords.shape must match.Found: %s and %s.'
        % (decoy_coords.shape, gt_coords.shape)
    )
  # Include and align all residues unless specified otherwise.
  if include_idxs is None:
    include_idxs = np.arange(decoy_coords.shape[0])
  if align_idxs is None:
    align_idxs = include_idxs
  aligned_decoy_coords = align(
      x=decoy_coords,
      y=gt_coords,
      x_indices=align_idxs,
      y_indices=align_idxs,
  )
  deviations = np.linalg.norm(
      aligned_decoy_coords[include_idxs] - gt_coords[include_idxs], axis=1
  )
  return deviations


def rmsd_from_coords(
    decoy_coords: np.ndarray,
    gt_coords: np.ndarray,
    align_idxs: np.ndarray | None = None,
    include_idxs: np.ndarray | None = None,
) -> float:
  """Computes the *aligned* RMSD of two Mx3 np arrays of coordinates.

  Args:
    decoy_coords: [M, 3] np array of decoy atom coordinates.
    gt_coords: [M, 3] np array of gt atom coordinates.
    align_idxs: [M] np array of indices specifying coordinates to align on.
      Defaults to None, in which case all the include_idx (see after) are used.
    include_idxs: [M] np array of indices specifying coordinates to score.
      Defaults to None, in which case all indices are used for scoring.

  Returns:
    rmsd value of the aligned decoy and gt coordinates.
  """
  deviations = deviations_from_coords(
      decoy_coords, gt_coords, align_idxs, include_idxs
  )
  return np.sqrt(np.mean(np.square(deviations)))
