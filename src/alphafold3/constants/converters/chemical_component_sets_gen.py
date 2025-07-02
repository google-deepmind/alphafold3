# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Script for updating chemical_component_sets.py."""

from collections.abc import Mapping, Sequence
import pathlib
import pickle
import re
import sys

from alphafold3.common import resources
import tqdm
from alphafold3.common import safe_pickle


_CCD_PICKLE_FILE = resources.filename(
    'constants/converters/ccd.pickle'
)


def find_ions_and_glycans_in_ccd(
    ccd: Mapping[str, Mapping[str, Sequence[str]]],
) -> dict[str, frozenset[str]]:
  """Finds glycans and ions in all version of CCD."""
  glycans_linking = []
  glycans_other = []
  ions = []
  for name, comp in tqdm.tqdm(ccd.items(), disable=None):
    if name == 'UNX':
      continue  # Skip "unknown atom or ion".
    comp_type = comp['_chem_comp.type'][0].lower()
    # Glycans have the type 'saccharide'.
    if re.findall(r'\bsaccharide\b', comp_type):
      # Separate out linking glycans from others.
      if 'linking' in comp_type:
        glycans_linking.append(name)
      else:
        glycans_other.append(name)

    # Ions have the word 'ion' in their name.
    comp_name = comp['_chem_comp.name'][0].lower()
    if re.findall(r'\bion\b', comp_name):
      ions.append(name)
  result = dict(
      glycans_linking=frozenset(glycans_linking),
      glycans_other=frozenset(glycans_other),
      ions=frozenset(ions),
  )

  return result


def main(argv: Sequence[str]) -> None:
  if len(argv) != 2:
    raise ValueError(
        'Directory to write to must be specified as a command-line arguments.'
    )

  print(f'Loading {_CCD_PICKLE_FILE}', flush=True)
  with open(_CCD_PICKLE_FILE, 'rb') as f:
    ccd: Mapping[str, Mapping[str, Sequence[str]]] = safe_pickle.safe_load(f)
  output_path = pathlib.Path(argv[1])
  output_path.parent.mkdir(exist_ok=True)
  print('Finding ions and glycans', flush=True)
  result = find_ions_and_glycans_in_ccd(ccd)
  print(f'writing to {output_path}', flush=True)
  with output_path.open('wb') as f:
    pickle.dump(result, f)
  print('Done', flush=True)


if __name__ == '__main__':
  main(sys.argv)
