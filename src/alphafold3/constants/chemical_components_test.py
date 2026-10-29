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
from alphafold3.constants import chemical_components


# 4NZ renamed to LIG-DM, ZZR modified to check precedence works.
_USER_CCD = """data_LIG-DM
# 
_chem_comp.id                                    LIG-DM 
_chem_comp.name                                  4-nitrobenzenesulfonamide 
_chem_comp.type                                  NON-POLYMER 
_chem_comp.pdbx_type                             HETAIN 
_chem_comp.formula                               "C6 H6 N2 O4 S" 
_chem_comp.mon_nstd_parent_comp_id               ? 
_chem_comp.pdbx_synonyms                         ? 
_chem_comp.pdbx_formal_charge                    0 
_chem_comp.pdbx_initial_date                     2010-04-01 
_chem_comp.pdbx_modified_date                    2011-06-04 
_chem_comp.pdbx_ambiguous_flag                   N 
_chem_comp.pdbx_release_status                   REL 
_chem_comp.pdbx_replaced_by                      ? 
_chem_comp.pdbx_replaces                         ? 
_chem_comp.formula_weight                        202.188 
_chem_comp.one_letter_code                       ? 
_chem_comp.three_letter_code                     LIG-DM 
_chem_comp.pdbx_model_coordinates_details        ? 
_chem_comp.pdbx_model_coordinates_missing_flag   N 
_chem_comp.pdbx_ideal_coordinates_details        Corina 
_chem_comp.pdbx_ideal_coordinates_missing_flag   N 
_chem_comp.pdbx_model_coordinates_db_code        3M8T 
_chem_comp.pdbx_subcomponent_list                ? 
_chem_comp.pdbx_processing_site                  PDBJ 
# 
loop_
_chem_comp_atom.comp_id 
_chem_comp_atom.atom_id 
_chem_comp_atom.alt_atom_id 
_chem_comp_atom.type_symbol 
_chem_comp_atom.charge 
_chem_comp_atom.pdbx_align 
_chem_comp_atom.pdbx_aromatic_flag 
_chem_comp_atom.pdbx_leaving_atom_flag 
_chem_comp_atom.pdbx_stereo_config 
_chem_comp_atom.model_Cartn_x 
_chem_comp_atom.model_Cartn_y 
_chem_comp_atom.model_Cartn_z 
_chem_comp_atom.pdbx_model_Cartn_x_ideal 
_chem_comp_atom.pdbx_model_Cartn_y_ideal 
_chem_comp_atom.pdbx_model_Cartn_z_ideal 
_chem_comp_atom.pdbx_component_atom_id 
_chem_comp_atom.pdbx_component_comp_id 
_chem_comp_atom.pdbx_ordinal 
LIG-DM S    S    S 0  1 N N N 4.225 10.331 -1.014 -2.542 -0.000 0.095  S    LIG-DM 1  
LIG-DM C1   C1   C 0  1 Y N N 7.695 12.230 0.978  1.984  -0.000 0.007  C1   LIG-DM 2  
LIG-DM N1   N1   N 0  1 N N N 4.373 8.737  -0.914 -3.062 0.002  -1.476 N1   LIG-DM 3  
LIG-DM O1   O1   O 0  1 N N N 3.016 10.900 -0.432 -2.932 1.249  0.649  O1   LIG-DM 4  
LIG-DM C2   C2   C 0  1 Y N N 7.656 10.848 0.833  1.293  1.197  0.021  C2   LIG-DM 5  
LIG-DM O2   O2   O 0  1 N N N 4.400 10.752 -2.379 -2.932 -1.250 0.646  O2   LIG-DM 6  
LIG-DM C3   C3   C 0  1 Y N N 6.563 10.256 0.194  -0.089 1.197  0.049  C3   LIG-DM 7  
LIG-DM C4   C4   C 0  1 Y N N 5.506 11.024 -0.260 -0.780 -0.000 0.061  C4   LIG-DM 8  
LIG-DM C5   C5   C 0  1 Y N N 5.559 12.410 -0.128 -0.089 -1.197 0.047  C5   LIG-DM 9  
LIG-DM C6   C6   C 0  1 Y N N 6.658 13.017 0.476  1.293  -1.197 0.020  C6   LIG-DM 10 
LIG-DM N7   N7   N 1  1 N N N 8.743 12.872 1.645  3.464  0.000  -0.029 N7   LIG-DM 11 
LIG-DM O12  O12  O 0  1 N N N 9.548 12.220 2.273  4.073  1.055  -0.045 O12  LIG-DM 12 
LIG-DM O13  O13  O -1 1 N N N 8.803 14.101 1.620  4.073  -1.055 -0.041 O13  LIG-DM 13 
LIG-DM HN1  HN1  H 0  1 N N N 5.232 8.454  -1.341 -2.418 0.002  -2.201 HN1  LIG-DM 14 
LIG-DM HN1A HN1A H 0  0 N N N 4.373 8.462  0.047  -4.012 0.002  -1.673 HN1A LIG-DM 15 
LIG-DM H2   H2   H 0  1 N N N 8.463 10.237 1.211  1.833  2.132  0.012  H2   LIG-DM 16 
LIG-DM H3   H3   H 0  1 N N N 6.543 9.185  0.053  -0.629 2.132  0.061  H3   LIG-DM 17 
LIG-DM H5   H5   H 0  1 N N N 4.744 13.017 -0.496 -0.629 -2.132 0.056  H5   LIG-DM 18 
LIG-DM H6   H6   H 0  1 N N N 6.707 14.093 0.555  1.833  -2.132 0.009  H6   LIG-DM 19 
# 
loop_
_chem_comp_bond.comp_id 
_chem_comp_bond.atom_id_1 
_chem_comp_bond.atom_id_2 
_chem_comp_bond.value_order 
_chem_comp_bond.pdbx_aromatic_flag 
_chem_comp_bond.pdbx_stereo_config 
_chem_comp_bond.pdbx_ordinal 
LIG-DM S  N1   SING N N 1  
LIG-DM S  O1   DOUB N N 2  
LIG-DM S  O2   DOUB N N 3  
LIG-DM S  C4   SING N N 4  
LIG-DM C1 C2   DOUB Y N 5  
LIG-DM C1 C6   SING Y N 6  
LIG-DM C1 N7   SING N N 7  
LIG-DM C2 C3   SING Y N 8  
LIG-DM C3 C4   DOUB Y N 9  
LIG-DM C4 C5   SING Y N 10 
LIG-DM C5 C6   DOUB Y N 11 
LIG-DM N7 O12  DOUB N N 12 
LIG-DM N7 O13  SING N N 13 
LIG-DM N1 HN1  SING N N 14 
LIG-DM N1 HN1A SING N N 15 
LIG-DM C2 H2   SING N N 16 
LIG-DM C3 H3   SING N N 17 
LIG-DM C5 H5   SING N N 18 
LIG-DM C6 H6   SING N N 19 
# 
loop_
_pdbx_chem_comp_descriptor.comp_id 
_pdbx_chem_comp_descriptor.type 
_pdbx_chem_comp_descriptor.program 
_pdbx_chem_comp_descriptor.program_version 
_pdbx_chem_comp_descriptor.descriptor 
LIG-DM SMILES           ACDLabs              12.01 "O=S(=O)(N)c1ccc([N+]([O-])=O)cc1"                                       
LIG-DM SMILES_CANONICAL CACTVS               3.370 "N[S](=O)(=O)c1ccc(cc1)[N+]([O-])=O"                                     
LIG-DM SMILES           CACTVS               3.370 "N[S](=O)(=O)c1ccc(cc1)[N+]([O-])=O"                                     
LIG-DM SMILES_CANONICAL "OpenEye OEToolkits" 1.7.0 "c1cc(ccc1[N+](=O)[O-])S(=O)(=O)N"                                       
LIG-DM SMILES           "OpenEye OEToolkits" 1.7.0 "c1cc(ccc1[N+](=O)[O-])S(=O)(=O)N"                                       
LIG-DM InChI            InChI                1.03  "InChI=1S/C6H6N2O4S/c7-13(11,12)6-3-1-5(2-4-6)8(9)10/h1-4H,(H2,7,11,12)" 
LIG-DM InChIKey         InChI                1.03  QWKKYJLAUWFPDB-UHFFFAOYSA-N                                              
# 
loop_
_pdbx_chem_comp_identifier.comp_id 
_pdbx_chem_comp_identifier.type 
_pdbx_chem_comp_identifier.program 
_pdbx_chem_comp_identifier.program_version 
_pdbx_chem_comp_identifier.identifier 
LIG-DM "SYSTEMATIC NAME" ACDLabs              12.01 4-nitrobenzenesulfonamide 
LIG-DM "SYSTEMATIC NAME" "OpenEye OEToolkits" 1.7.0 4-nitrobenzenesulfonamide 
# 
loop_
_pdbx_chem_comp_audit.comp_id 
_pdbx_chem_comp_audit.action_type 
_pdbx_chem_comp_audit.date 
_pdbx_chem_comp_audit.processing_site 
LIG-DM "Create component"     2010-04-01 PDBJ 
LIG-DM "Modify aromatic_flag" 2011-06-04 RCSB 
LIG-DM "Modify descriptor"    2011-06-04 RCSB 
# 

data_ZZR
# 
_chem_comp.id                                    ZZR 
_chem_comp.name                                  "Test name" 
_chem_comp.type                                  non-polymer 
_chem_comp.pdbx_type                             HETAIN 
_chem_comp.formula                               "C3 H6 N7" 
_chem_comp.mon_nstd_parent_comp_id               ? 
_chem_comp.pdbx_synonyms                         ? 
_chem_comp.pdbx_formal_charge                    1 
_chem_comp.pdbx_initial_date                     2009-11-10 
_chem_comp.pdbx_modified_date                    2011-06-04 
_chem_comp.pdbx_ambiguous_flag                   N 
_chem_comp.pdbx_release_status                   REL 
_chem_comp.pdbx_replaced_by                      ? 
_chem_comp.pdbx_replaces                         ? 
_chem_comp.formula_weight                        140.127 
_chem_comp.one_letter_code                       ? 
_chem_comp.three_letter_code                     ZZR 
_chem_comp.pdbx_model_coordinates_details        ? 
_chem_comp.pdbx_model_coordinates_missing_flag   N 
_chem_comp.pdbx_ideal_coordinates_details        Corina 
_chem_comp.pdbx_ideal_coordinates_missing_flag   N 
_chem_comp.pdbx_model_coordinates_db_code        2WXS 
_chem_comp.pdbx_subcomponent_list                ? 
_chem_comp.pdbx_processing_site                  EBI 
# 
loop_
_chem_comp_atom.comp_id 
_chem_comp_atom.atom_id 
_chem_comp_atom.alt_atom_id 
_chem_comp_atom.type_symbol 
_chem_comp_atom.charge 
_chem_comp_atom.pdbx_align 
_chem_comp_atom.pdbx_aromatic_flag 
_chem_comp_atom.pdbx_leaving_atom_flag 
_chem_comp_atom.pdbx_stereo_config 
_chem_comp_atom.model_Cartn_x 
_chem_comp_atom.model_Cartn_y 
_chem_comp_atom.model_Cartn_z 
_chem_comp_atom.pdbx_model_Cartn_x_ideal 
_chem_comp_atom.pdbx_model_Cartn_y_ideal 
_chem_comp_atom.pdbx_model_Cartn_z_ideal 
_chem_comp_atom.pdbx_component_atom_id 
_chem_comp_atom.pdbx_component_comp_id 
_chem_comp_atom.pdbx_ordinal 
ZZR NAC  NAC  N 0 1 Y N N -11.751 -2.583 -6.221  2.175  -0.751 -0.000 NAC  ZZR 1  
ZZR NAF  NAF  N 0 1 Y N N -12.220 -3.054 -7.283  1.126  -1.683 -0.000 NAF  ZZR 2  
ZZR CAI  CAI  C 0 1 Y N N -12.066 -2.158 -8.268  -0.035 -0.955 0.000  CAI  ZZR 3  
ZZR NAD  NAD  N 0 1 Y N N -12.289 -1.938 -9.566  -1.374 -1.053 -0.000 NAD  ZZR 4  
ZZR CAG  CAG  C 0 1 Y N N -11.873 -0.694 -9.770  -1.908 0.148  -0.000 CAG  ZZR 5  
ZZR NAA  NAA  N 0 1 N N N -11.934 -0.083 -10.956 -3.259 0.434  0.000  NAA  ZZR 6  
ZZR NAE  NAE  N 0 1 Y N N -11.363 -0.185 -8.623  -0.900 1.060  -0.000 NAE  ZZR 7  
ZZR NAJ  NAJ  N 1 1 Y N N -11.496 -1.054 -7.737  0.299  0.328  0.000  NAJ  ZZR 8  
ZZR CAH  CAH  C 0 1 Y N N -11.297 -1.349 -6.456  1.636  0.442  0.000  CAH  ZZR 9  
ZZR NAB  NAB  N 0 1 N N N -10.780 -0.537 -5.534  2.344  1.647  -0.000 NAB  ZZR 10 
ZZR HAF  HAF  H 0 1 N N N -12.642 -3.955 -7.380  1.207  -2.650 -0.000 HAF  ZZR 11 
ZZR HAA1 HAA1 H 0 0 N N N -11.949 0.908  -10.821 -3.907 -0.287 0.000  HAA1 ZZR 12 
ZZR HAA2 HAA2 H 0 0 N N N -11.134 -0.333 -11.501 -3.559 1.357  0.000  HAA2 ZZR 13 
ZZR HAE  HAE  H 0 1 N N N -10.955 0.721  -8.507  -0.990 2.026  -0.001 HAE  ZZR 14 
ZZR HAB1 HAB1 H 0 0 N N N -10.768 -1.011 -4.653  3.314  1.641  -0.001 HAB1 ZZR 15 
ZZR HAB2 HAB2 H 0 0 N N N -10.456 0.395  -5.697  1.865  2.491  0.000  HAB2 ZZR 16 
# 
loop_
_chem_comp_bond.comp_id 
_chem_comp_bond.atom_id_1 
_chem_comp_bond.atom_id_2 
_chem_comp_bond.value_order 
_chem_comp_bond.pdbx_aromatic_flag 
_chem_comp_bond.pdbx_stereo_config 
_chem_comp_bond.pdbx_ordinal 
ZZR NAC NAF  SING Y N 1  
ZZR NAC CAH  DOUB Y N 2  
ZZR NAF CAI  SING Y N 3  
ZZR CAI NAD  SING Y N 4  
ZZR CAI NAJ  DOUB Y N 5  
ZZR NAD CAG  DOUB Y N 6  
ZZR CAG NAA  SING N N 7  
ZZR CAG NAE  SING Y N 8  
ZZR NAE NAJ  SING Y N 9  
ZZR NAJ CAH  SING Y N 10 
ZZR CAH NAB  SING N N 11 
ZZR NAF HAF  SING N N 12 
ZZR NAA HAA1 SING N N 13 
ZZR NAA HAA2 SING N N 14 
ZZR NAE HAE  SING N N 15 
ZZR NAB HAB1 SING N N 16 
ZZR NAB HAB2 SING N N 17 
# 
loop_
_pdbx_chem_comp_descriptor.comp_id 
_pdbx_chem_comp_descriptor.type 
_pdbx_chem_comp_descriptor.program 
_pdbx_chem_comp_descriptor.program_version 
_pdbx_chem_comp_descriptor.descriptor 
ZZR SMILES           ACDLabs              10.04 "n1c(N)n[n+]2c(nnc12)N"                                           
ZZR SMILES_CANONICAL CACTVS               3.352 "Nc1[nH][n+]2c(N)n[nH]c2n1"                                       
ZZR SMILES           CACTVS               3.352 "Nc1[nH][n+]2c(N)n[nH]c2n1"                                       
ZZR SMILES_CANONICAL "OpenEye OEToolkits" 1.6.1 "c1([nH][n+]2c(n[nH]c2n1)N)N"                                     
ZZR SMILES           "OpenEye OEToolkits" 1.6.1 "c1([nH][n+]2c(n[nH]c2n1)N)N"                                     
ZZR InChI            InChI                1.03  "InChI=1S/C3H5N7/c4-1-6-3-8-7-2(5)10(3)9-1/h(H5,4,5,6,7,8,9)/p+1" 
ZZR InChIKey         InChI                1.03  JFMHNXXJHXPFLE-UHFFFAOYSA-O                                       
# 
loop_
_pdbx_chem_comp_identifier.comp_id 
_pdbx_chem_comp_identifier.type 
_pdbx_chem_comp_identifier.program 
_pdbx_chem_comp_identifier.program_version 
_pdbx_chem_comp_identifier.identifier 
ZZR "SYSTEMATIC NAME" ACDLabs              10.04 "3,6-diamino-1,5-dihydro[1,2,4]triazolo[4,3-b][1,2,4]triazol-4-ium"  
ZZR "SYSTEMATIC NAME" "OpenEye OEToolkits" 1.6.1 "1,5-dihydro-[1,2,4]triazolo[5,4-e][1,2,4]triazol-4-ium-3,6-diamine" 
# 
loop_
_pdbx_chem_comp_audit.comp_id 
_pdbx_chem_comp_audit.action_type 
_pdbx_chem_comp_audit.date 
_pdbx_chem_comp_audit.processing_site 
ZZR "Create component"     2009-11-10 EBI  
ZZR "Modify aromatic_flag" 2011-06-04 RCSB 
ZZR "Modify descriptor"    2011-06-04 RCSB 
#  
"""


class ChemicalComponentDatabaseTest(parameterized.TestCase):

  def test_user_ccd(self):
    ccd = chemical_components.Ccd(ccd_msgpack_path=None, user_ccd=_USER_CCD)
    # 4NZ renamed to LIG-DM, ZZR modified to check precedence works.
    self.assertIsNotNone(ccd.get('LIG-DM'))
    entry_4nz = dict(ccd['4NZ'])
    entry_lig_dm = dict(ccd['LIG-DM'])
    for key in (
        '_chem_comp.id',
        '_chem_comp.three_letter_code',
        '_chem_comp_atom.comp_id',
        '_chem_comp_atom.pdbx_component_comp_id',
        '_chem_comp_bond.comp_id',
        '_pdbx_chem_comp_audit.comp_id',
        '_pdbx_chem_comp_descriptor.comp_id',
        '_pdbx_chem_comp_identifier.comp_id',
        'data_',
    ):
      del entry_4nz[key]
      del entry_lig_dm[key]
    self.assertEqual(entry_4nz, entry_lig_dm)

    zzr_comp = ccd.get('ZZR')
    self.assertIsNotNone(zzr_comp)
    self.assertEqual(zzr_comp['_chem_comp.name'][0], 'Test name')

  def test_bad_user_ccd(self):
    with self.assertRaises(ValueError):
      chemical_components.Ccd(ccd_msgpack_path=None, user_ccd='bla bla')

  def test_arg(self):
    ccd = chemical_components.Ccd()
    self.assertSequenceEqual(ccd['ARG']['data_'], ['ARG'])
    self.assertSequenceEqual(
        ccd['ARG']['_chem_comp.type'], ['L-PEPTIDE LINKING']
    )
    self.assertSequenceEqual(
        ccd['ARG']['_pdbx_chem_comp_descriptor.descriptor'],
        [
            'O=C(O)C(N)CCCN\\C(=[NH2+])N',
            'N[C@@H](CCCNC(N)=[NH2+])C(O)=O',
            'N[CH](CCCNC(N)=[NH2+])C(O)=O',
            'C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N',
            'C(CC(C(=O)O)N)CNC(=[NH2+])N',
            'InChI=1S/C6H14N4O2/c7-4(5(11)12)2-1-3-10-6(8)9/h4H,1-3,7H2,(H,11,12)(H4,8,9,10)/p+1/t4-/m0/s1',
            'ODKSFYDXXFIFQN-BYPYZUCNSA-O',
        ],
    )

  def test_arg_smiles(self):
    arg_ccd = chemical_components.Ccd()['ARG']
    component_info = chemical_components.mmcif_to_info(arg_ccd)
    self.assertEqual(
        component_info.pdbx_smiles, 'C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N'
    )

  @parameterized.parameters(
      dict(res_name='ARG', atom_name='CA', exp_type_symbol='C'),
      dict(res_name='FOO', atom_name='BAR', exp_type_symbol='?'),
      dict(res_name='ARG', atom_name='BAR', exp_type_symbol='?'),
  )
  def test_type_symbol(self, res_name, atom_name, exp_type_symbol):
    ccd = chemical_components.Ccd()
    self.assertEqual(
        chemical_components.type_symbol(
            ccd, res_name=res_name, atom_name=atom_name
        ),
        exp_type_symbol,
    )


if __name__ == '__main__':
  absltest.main()
