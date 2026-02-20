# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Script for building intermediate data."""

from importlib import resources
import pathlib
import site
import sysconfig

import alphafold3.constants.converters
from alphafold3.constants.converters import ccd_pickle_gen
from alphafold3.constants.converters import chemical_component_sets_gen


def build_data():
  """Builds intermediate data."""
  for site_path in site.getsitepackages():
    path = pathlib.Path(site_path) / 'share/libcifpp/components.cif'
    if path.exists():
      cif_path = path
      break
  else:
    # Fallback: check sysconfig data path for installations outside Python prefix.
    data_path = pathlib.Path(sysconfig.get_path('data')) / 'share/libcifpp/components.cif'
    if data_path.exists():
      cif_path = data_path
    else:
      import alphafold3
      af3_data_path = pathlib.Path(alphafold3.__path__[0]) / '../share/libcifpp/components.cif'
      if af3_data_path.exists():
        cif_path = af3_data_path
      else:
        raise ValueError('Could not find components.cif')

  out_root = resources.files(alphafold3.constants.converters)
  ccd_pickle_path = out_root.joinpath('ccd.pickle')
  chemical_component_sets_pickle_path = out_root.joinpath(
      'chemical_component_sets.pickle'
  )
  ccd_pickle_gen.main(['', str(cif_path), str(ccd_pickle_path)])
  chemical_component_sets_gen.main(
      ['', str(chemical_component_sets_pickle_path)]
  )
