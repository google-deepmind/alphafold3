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

"""Constants associated with residue names."""

from collections.abc import Mapping
import functools
import sys

# pyformat: disable
# common_typos_disable
CCD_NAME_TO_ONE_LETTER: Mapping[str, str] = {
    '00C': 'C', '01W': 'X', '02K': 'A', '03Y': 'C', '07O': 'C', '08P': 'C',
    '0A0': 'D', '0A1': 'Y', '0A2': 'K', '0A8': 'C', '0AA': 'V', '0AB': 'V',
    '0AC': 'G', '0AD': 'G', '0AF': 'W', '0AG': 'L', '0AH': 'S', '0AK': 'D',
    '0AM': 'A', '0AP': 'C', '0AU': 'U', '0AV': 'A', '0AZ': 'P', '0BN': 'F',
    '0C': 'C', '0CS': 'A', '0DC': 'C', '0DG': 'G', '0DT': 'T', '0FL': 'A',
    '0G': 'G', '0NC': 'A', '0SP': 'A', '0U': 'U', '10C': 'C', '125': 'U',
    '126': 'U', '127': 'U', '128': 'N', '12A': 'A', '143': 'C', '193': 'X',
    '1AP': 'A', '1MA': 'A', '1MG': 'G', '1PA': 'F', '1PI': 'A', '1PR': 'N',
    '1SC': 'C', '1TQ': 'W', '1TY': 'Y', '1X6': 'S', '200': 'F', '23F': 'F',
    '23S': 'X', '26B': 'T', '2AD': 'X', '2AG': 'A', '2AO': 'X', '2AR': 'A',
    '2AS': 'X', '2AT': 'T', '2AU': 'U', '2BD': 'I', '2BT': 'T', '2BU': 'A',
    '2CO': 'C', '2DA': 'A', '2DF': 'N', '2DM': 'N', '2DO': 'X', '2DT': 'T',
    '2EG': 'G', '2FE': 'N', '2FI': 'N', '2FM': 'M', '2GT': 'T', '2HF': 'H',
    '2LU': 'L', '2MA': 'A', '2MG': 'G', '2ML': 'L', '2MR': 'R', '2MT': 'P',
    '2MU': 'U', '2NT': 'T', '2OM': 'U', '2OT': 'T', '2PI': 'X', '2PR': 'G',
    '2SA': 'N', '2SI': 'X', '2ST': 'T', '2TL': 'T', '2TY': 'Y', '2VA': 'V',
    '2XA': 'C', '32S': 'X', '32T': 'X', '3AH': 'H', '3AR': 'X', '3CF': 'F',
    '3DA': 'A', '3DR': 'N', '3GA': 'A', '3MD': 'D', '3ME': 'U', '3NF': 'Y',
    '3QN': 'K', '3TY': 'X', '3XH': 'G', '4AC': 'N', '4BF': 'Y', '4CF': 'F',
    '4CY': 'M', '4DP': 'W', '4FB': 'P', '4FW': 'W', '4HT': 'W', '4IN': 'W',
    '4MF': 'N', '4MM': 'X', '4OC': 'C', '4PC': 'C', '4PD': 'C', '4PE': 'C',
    '4PH': 'F', '4SC': 'C', '4SU': 'U', '4TA': 'N', '4U7': 'A', '56A': 'H',
    '5AA': 'A', '5AB': 'A', '5AT': 'T', '5BU': 'U', '5CG': 'G', '5CM': 'C',
    '5CS': 'C', '5FA': 'A', '5FC': 'C', '5FU': 'U', '5HP': 'E', '5HT': 'T',
    '5HU': 'U', '5IC': 'C', '5IT': 'T', '5IU': 'U', '5MC': 'C', '5MD': 'N',
    '5MU': 'U', '5NC': 'C', '5PC': 'C', '5PY': 'T', '5SE': 'U', '64T': 'T',
    '6CL': 'K', '6CT': 'T', '6CW': 'W', '6HA': 'A', '6HC': 'C', '6HG': 'G',
    '6HN': 'K', '6HT': 'T', '6IA': 'A', '6MA': 'A', '6MC': 'A', '6MI': 'N',
    '6MT': 'A', '6MZ': 'N', '6OG': 'G', '70U': 'U', '7DA': 'A', '7GU': 'G',
    '7JA': 'I', '7MG': 'G', '8AN': 'A', '8FG': 'G', '8MG': 'G', '8OG': 'G',
    '9NE': 'E', '9NF': 'F', '9NR': 'R', '9NV': 'V', 'A': 'A', 'A1P': 'N',
    'A23': 'A', 'A2L': 'A', 'A2M': 'A', 'A34': 'A', 'A35': 'A', 'A38': 'A',
    'A39': 'A', 'A3A': 'A', 'A3P': 'A', 'A40': 'A', 'A43': 'A', 'A44': 'A',
    'A47': 'A', 'A5L': 'A', 'A5M': 'C', 'A5N': 'N', 'A5O': 'A', 'A66': 'X',
    'AA3': 'A', 'AA4': 'A', 'AAR': 'R', 'AB7': 'X', 'ABA': 'A', 'ABR': 'A',
    'ABS': 'A', 'ABT': 'N', 'ACB': 'D', 'ACL': 'R', 'AD2': 'A', 'ADD': 'X',
    'ADX': 'N', 'AEA': 'X', 'AEI': 'D', 'AET': 'A', 'AFA': 'N', 'AFF': 'N',
    'AFG': 'G', 'AGM': 'R', 'AGT': 'C', 'AHB': 'N', 'AHH': 'X', 'AHO': 'A',
    'AHP': 'A', 'AHS': 'X', 'AHT': 'X', 'AIB': 'A', 'AKL': 'D', 'AKZ': 'D',
    'ALA': 'A', 'ALC': 'A', 'ALM': 'A', 'ALN': 'A', 'ALO': 'T', 'ALQ': 'X',
    'ALS': 'A', 'ALT': 'A', 'ALV': 'A', 'ALY': 'K', 'AN8': 'A', 'AP7': 'A',
    'APE': 'X', 'APH': 'A', 'API': 'K', 'APK': 'K', 'APM': 'X', 'APP': 'X',
    'AR2': 'R', 'AR4': 'E', 'AR7': 'R', 'ARG': 'R', 'ARM': 'R', 'ARO': 'R',
    'ARV': 'X', 'AS': 'A', 'AS2': 'D', 'AS9': 'X', 'ASA': 'D', 'ASB': 'D',
    'ASI': 'D', 'ASK': 'D', 'ASL': 'D', 'ASM': 'X', 'ASN': 'N', 'ASP': 'D',
    'ASQ': 'D', 'ASU': 'N', 'ASX': 'B', 'ATD': 'T', 'ATL': 'T', 'ATM': 'T',
    'AVC': 'A', 'AVN': 'X', 'AYA': 'A', 'AZK': 'K', 'AZS': 'S', 'AZY': 'Y',
    'B1F': 'F', 'B1P': 'N', 'B2A': 'A', 'B2F': 'F', 'B2I': 'I', 'B2V': 'V',
    'B3A': 'A', 'B3D': 'D', 'B3E': 'E', 'B3K': 'K', 'B3L': 'X', 'B3M': 'X',
    'B3Q': 'X', 'B3S': 'S', 'B3T': 'X', 'B3U': 'H', 'B3X': 'N', 'B3Y': 'Y',
    'BB6': 'C', 'BB7': 'C', 'BB8': 'F', 'BB9': 'C', 'BBC': 'C', 'BCS': 'C',
    'BE2': 'X', 'BFD': 'D', 'BG1': 'S', 'BGM': 'G', 'BH2': 'D', 'BHD': 'D',
    'BIF': 'F', 'BIL': 'X', 'BIU': 'I', 'BJH': 'X', 'BLE': 'L', 'BLY': 'K',
    'BMP': 'N', 'BMT': 'T', 'BNN': 'F', 'BNO': 'X', 'BOE': 'T', 'BOR': 'R',
    'BPE': 'C', 'BRU': 'U', 'BSE': 'S', 'BT5': 'N', 'BTA': 'L', 'BTC': 'C',
    'BTR': 'W', 'BUC': 'C', 'BUG': 'V', 'BVP': 'U', 'BZG': 'N', 'C': 'C',
    'C1X': 'K', 'C25': 'C', 'C2L': 'C', 'C2S': 'C', 'C31': 'C', 'C32': 'C',
    'C34': 'C', 'C36': 'C', 'C37': 'C', 'C38': 'C', 'C3Y': 'C', 'C42': 'C',
    'C43': 'C', 'C45': 'C', 'C46': 'C', 'C49': 'C', 'C4R': 'C', 'C4S': 'C',
    'C5C': 'C', 'C66': 'X', 'C6C': 'C', 'CAF': 'C', 'CAL': 'X', 'CAR': 'C',
    'CAS': 'C', 'CAV': 'X', 'CAY': 'C', 'CB2': 'C', 'CBR': 'C', 'CBV': 'C',
    'CCC': 'C', 'CCL': 'K', 'CCS': 'C', 'CDE': 'X', 'CDV': 'X', 'CDW': 'C',
    'CEA': 'C', 'CFL': 'C', 'CG1': 'G', 'CGA': 'E', 'CGU': 'E', 'CH': 'C',
    'CHF': 'X', 'CHG': 'X', 'CHP': 'G', 'CHS': 'X', 'CIR': 'R', 'CLE': 'L',
    'CLG': 'K', 'CLH': 'K', 'CM0': 'N', 'CME': 'C', 'CMH': 'C', 'CML': 'C',
    'CMR': 'C', 'CMT': 'C', 'CNU': 'U', 'CP1': 'C', 'CPC': 'X', 'CPI': 'X',
    'CR5': 'G', 'CS0': 'C', 'CS1': 'C', 'CS3': 'C', 'CS4': 'C', 'CS8': 'N',
    'CSA': 'C', 'CSB': 'C', 'CSD': 'C', 'CSE': 'C', 'CSF': 'C', 'CSI': 'G',
    'CSJ': 'C', 'CSL': 'C', 'CSO': 'C', 'CSP': 'C', 'CSR': 'C', 'CSS': 'C',
    'CSU': 'C', 'CSW': 'C', 'CSX': 'C', 'CSZ': 'C', 'CTE': 'W', 'CTG': 'T',
    'CTH': 'T', 'CUC': 'X', 'CWR': 'S', 'CXM': 'M', 'CY0': 'C', 'CY1': 'C',
    'CY3': 'C', 'CY4': 'C', 'CYA': 'C', 'CYD': 'C', 'CYF': 'C', 'CYG': 'C',
    'CYJ': 'X', 'CYM': 'C', 'CYQ': 'C', 'CYR': 'C', 'CYS': 'C', 'CZ2': 'C',
    'CZZ': 'C', 'D11': 'T', 'D1P': 'N', 'D3': 'N', 'D33': 'N', 'D3P': 'G',
    'D3T': 'T', 'D4M': 'T', 'D4P': 'X', 'DA': 'A', 'DA2': 'X', 'DAB': 'A',
    'DAH': 'F', 'DAL': 'A', 'DAR': 'R', 'DAS': 'D', 'DBB': 'T', 'DBM': 'N',
    'DBS': 'S', 'DBU': 'T', 'DBY': 'Y', 'DBZ': 'A', 'DC': 'C', 'DC2': 'C',
    'DCG': 'G', 'DCI': 'X', 'DCL': 'X', 'DCT': 'C', 'DCY': 'C', 'DDE': 'H',
    'DDG': 'G', 'DDN': 'U', 'DDX': 'N', 'DFC': 'C', 'DFG': 'G', 'DFI': 'X',
    'DFO': 'X', 'DFT': 'N', 'DG': 'G', 'DGH': 'G', 'DGI': 'G', 'DGL': 'E',
    'DGN': 'Q', 'DHA': 'S', 'DHI': 'H', 'DHL': 'X', 'DHN': 'V', 'DHP': 'X',
    'DHU': 'U', 'DHV': 'V', 'DI': 'I', 'DIL': 'I', 'DIR': 'R', 'DIV': 'V',
    'DLE': 'L', 'DLS': 'K', 'DLY': 'K', 'DM0': 'K', 'DMH': 'N', 'DMK': 'D',
    'DMT': 'X', 'DN': 'N', 'DNE': 'L', 'DNG': 'L', 'DNL': 'K', 'DNM': 'L',
    'DNP': 'A', 'DNR': 'C', 'DNS': 'K', 'DOA': 'X', 'DOC': 'C', 'DOH': 'D',
    'DON': 'L', 'DPB': 'T', 'DPH': 'F', 'DPL': 'P', 'DPP': 'A', 'DPQ': 'Y',
    'DPR': 'P', 'DPY': 'N', 'DRM': 'U', 'DRP': 'N', 'DRT': 'T', 'DRZ': 'N',
    'DSE': 'S', 'DSG': 'N', 'DSN': 'S', 'DSP': 'D', 'DT': 'T', 'DTH': 'T',
    'DTR': 'W', 'DTY': 'Y', 'DU': 'U', 'DVA': 'V', 'DXD': 'N', 'DXN': 'N',
    'DYS': 'C', 'DZM': 'A', 'E': 'A', 'E1X': 'A', 'ECC': 'Q', 'EDA': 'A',
    'EFC': 'C', 'EHP': 'F', 'EIT': 'T', 'ENP': 'N', 'ESB': 'Y', 'ESC': 'M',
    'EXB': 'X', 'EXY': 'L', 'EY5': 'N', 'EYS': 'X', 'F2F': 'F', 'FA2': 'A',
    'FA5': 'N', 'FAG': 'N', 'FAI': 'N', 'FB5': 'A', 'FB6': 'A', 'FCL': 'F',
    'FFD': 'N', 'FGA': 'E', 'FGL': 'G', 'FGP': 'S', 'FHL': 'X', 'FHO': 'K',
    'FHU': 'U', 'FLA': 'A', 'FLE': 'L', 'FLT': 'Y', 'FME': 'M', 'FMG': 'G',
    'FMU': 'N', 'FOE': 'C', 'FOX': 'G', 'FP9': 'P', 'FPA': 'F', 'FRD': 'X',
    'FT6': 'W', 'FTR': 'W', 'FTY': 'Y', 'FVA': 'V', 'FZN': 'K', 'G': 'G',
    'G25': 'G', 'G2L': 'G', 'G2S': 'G', 'G31': 'G', 'G32': 'G', 'G33': 'G',
    'G36': 'G', 'G38': 'G', 'G42': 'G', 'G46': 'G', 'G47': 'G', 'G48': 'G',
    'G49': 'G', 'G4P': 'N', 'G7M': 'G', 'GAO': 'G', 'GAU': 'E', 'GCK': 'C',
    'GCM': 'X', 'GDP': 'G', 'GDR': 'G', 'GFL': 'G', 'GGL': 'E', 'GH3': 'G',
    'GHG': 'Q', 'GHP': 'G', 'GL3': 'G', 'GLH': 'Q', 'GLJ': 'E', 'GLK': 'E',
    'GLM': 'X', 'GLN': 'Q', 'GLQ': 'E', 'GLU': 'E', 'GLX': 'Z', 'GLY': 'G',
    'GLZ': 'G', 'GMA': 'E', 'GMS': 'G', 'GMU': 'U', 'GN7': 'G', 'GND': 'X',
    'GNE': 'N', 'GOM': 'G', 'GPL': 'K', 'GS': 'G', 'GSC': 'G', 'GSR': 'G',
    'GSS': 'G', 'GSU': 'E', 'GT9': 'C', 'GTP': 'G', 'GVL': 'X', 'H2U': 'U',
    'H5M': 'P', 'HAC': 'A', 'HAR': 'R', 'HBN': 'H', 'HCS': 'X', 'HDP': 'U',
    'HEU': 'U', 'HFA': 'X', 'HGL': 'X', 'HHI': 'H', 'HIA': 'H', 'HIC': 'H',
    'HIP': 'H', 'HIQ': 'H', 'HIS': 'H', 'HL2': 'L', 'HLU': 'L', 'HMR': 'R',
    'HOL': 'N', 'HPC': 'F', 'HPE': 'F', 'HPH': 'F', 'HPQ': 'F', 'HQA': 'A',
    'HRG': 'R', 'HRP': 'W', 'HS8': 'H', 'HS9': 'H', 'HSE': 'S', 'HSL': 'S',
    'HSO': 'H', 'HTI': 'C', 'HTN': 'N', 'HTR': 'W', 'HV5': 'A', 'HVA': 'V',
    'HY3': 'P', 'HYP': 'P', 'HZP': 'P', 'I': 'I', 'I2M': 'I', 'I58': 'K',
    'I5C': 'C', 'IAM': 'A', 'IAR': 'R', 'IAS': 'D', 'IC': 'C', 'IEL': 'K',
    'IG': 'G', 'IGL': 'G', 'IGU': 'G', 'IIL': 'I', 'ILE': 'I', 'ILG': 'E',
    'ILX': 'I', 'IMC': 'C', 'IML': 'I', 'IOY': 'F', 'IPG': 'G', 'IPN': 'N',
    'IRN': 'N', 'IT1': 'K', 'IU': 'U', 'IYR': 'Y', 'IYT': 'T', 'IZO': 'M',
    'JJJ': 'C', 'JJK': 'C', 'JJL': 'C', 'JW5': 'N', 'K1R': 'C', 'KAG': 'G',
    'KCX': 'K', 'KGC': 'K', 'KNB': 'A', 'KOR': 'M', 'KPI': 'K', 'KST': 'K',
    'KYQ': 'K', 'L2A': 'X', 'LA2': 'K', 'LAA': 'D', 'LAL': 'A', 'LBY': 'K',
    'LC': 'C', 'LCA': 'A', 'LCC': 'N', 'LCG': 'G', 'LCH': 'N', 'LCK': 'K',
    'LCX': 'K', 'LDH': 'K', 'LED': 'L', 'LEF': 'L', 'LEH': 'L', 'LEI': 'V',
    'LEM': 'L', 'LEN': 'L', 'LET': 'X', 'LEU': 'L', 'LEX': 'L', 'LG': 'G',
    'LGP': 'G', 'LHC': 'X', 'LHU': 'U', 'LKC': 'N', 'LLP': 'K', 'LLY': 'K',
    'LME': 'E', 'LMF': 'K', 'LMQ': 'Q', 'LMS': 'N', 'LP6': 'K', 'LPD': 'P',
    'LPG': 'G', 'LPL': 'X', 'LPS': 'S', 'LSO': 'X', 'LTA': 'X', 'LTR': 'W',
    'LVG': 'G', 'LVN': 'V', 'LYF': 'K', 'LYK': 'K', 'LYM': 'K', 'LYN': 'K',
    'LYR': 'K', 'LYS': 'K', 'LYX': 'K', 'LYZ': 'K', 'M0H': 'C', 'M1G': 'G',
    'M2G': 'G', 'M2L': 'K', 'M2S': 'M', 'M30': 'G', 'M3L': 'K', 'M5M': 'C',
    'MA': 'A', 'MA6': 'A', 'MA7': 'A', 'MAA': 'A', 'MAD': 'A', 'MAI': 'R',
    'MBQ': 'Y', 'MBZ': 'N', 'MC1': 'S', 'MCG': 'X', 'MCL': 'K', 'MCS': 'C',
    'MCY': 'C', 'MD3': 'C', 'MD6': 'G', 'MDH': 'X', 'MDR': 'N', 'MEA': 'F',
    'MED': 'M', 'MEG': 'E', 'MEN': 'N', 'MEP': 'U', 'MEQ': 'Q', 'MET': 'M',
    'MEU': 'G', 'MF3': 'X', 'MG1': 'G', 'MGG': 'R', 'MGN': 'Q', 'MGQ': 'A',
    'MGV': 'G', 'MGY': 'G', 'MHL': 'L', 'MHO': 'M', 'MHS': 'H', 'MIA': 'A',
    'MIS': 'S', 'MK8': 'L', 'ML3': 'K', 'MLE': 'L', 'MLL': 'L', 'MLY': 'K',
    'MLZ': 'K', 'MME': 'M', 'MMO': 'R', 'MMT': 'T', 'MND': 'N', 'MNL': 'L',
    'MNU': 'U', 'MNV': 'V', 'MOD': 'X', 'MP8': 'P', 'MPH': 'X', 'MPJ': 'X',
    'MPQ': 'G', 'MRG': 'G', 'MSA': 'G', 'MSE': 'M', 'MSL': 'M', 'MSO': 'M',
    'MSP': 'X', 'MT2': 'M', 'MTR': 'T', 'MTU': 'A', 'MTY': 'Y', 'MVA': 'V',
    'N': 'N', 'N10': 'S', 'N2C': 'X', 'N5I': 'N', 'N5M': 'C', 'N6G': 'G',
    'N7P': 'P', 'NA8': 'A', 'NAL': 'A', 'NAM': 'A', 'NB8': 'N', 'NBQ': 'Y',
    'NC1': 'S', 'NCB': 'A', 'NCX': 'N', 'NCY': 'X', 'NDF': 'F', 'NDN': 'U',
    'NEM': 'H', 'NEP': 'H', 'NF2': 'N', 'NFA': 'F', 'NHL': 'E', 'NIT': 'X',
    'NIY': 'Y', 'NLE': 'L', 'NLN': 'L', 'NLO': 'L', 'NLP': 'L', 'NLQ': 'Q',
    'NMC': 'G', 'NMM': 'R', 'NMS': 'T', 'NMT': 'T', 'NNH': 'R', 'NP3': 'N',
    'NPH': 'C', 'NPI': 'A', 'NSK': 'X', 'NTY': 'Y', 'NVA': 'V', 'NYM': 'N',
    'NYS': 'C', 'NZH': 'H', 'O12': 'X', 'O2C': 'N', 'O2G': 'G', 'OAD': 'N',
    'OAS': 'S', 'OBF': 'X', 'OBS': 'X', 'OCS': 'C', 'OCY': 'C', 'ODP': 'N',
    'OHI': 'H', 'OHS': 'D', 'OIC': 'X', 'OIP': 'I', 'OLE': 'X', 'OLT': 'T',
    'OLZ': 'S', 'OMC': 'C', 'OMG': 'G', 'OMT': 'M', 'OMU': 'U', 'ONE': 'U',
    'ONH': 'A', 'ONL': 'X', 'OPR': 'R', 'ORN': 'A', 'ORQ': 'R', 'OSE': 'S',
    'OTB': 'X', 'OTH': 'T', 'OTY': 'Y', 'OXX': 'D', 'P': 'G', 'P1L': 'C',
    'P1P': 'N', 'P2T': 'T', 'P2U': 'U', 'P2Y': 'P', 'P5P': 'A', 'PAQ': 'Y',
    'PAS': 'D', 'PAT': 'W', 'PAU': 'A', 'PBB': 'C', 'PBF': 'F', 'PBT': 'N',
    'PCA': 'E', 'PCC': 'P', 'PCE': 'X', 'PCS': 'F', 'PDL': 'X', 'PDU': 'U',
    'PEC': 'C', 'PF5': 'F', 'PFF': 'F', 'PFX': 'X', 'PG1': 'S', 'PG7': 'G',
    'PG9': 'G', 'PGL': 'X', 'PGN': 'G', 'PGP': 'G', 'PGY': 'G', 'PHA': 'F',
    'PHD': 'D', 'PHE': 'F', 'PHI': 'F', 'PHL': 'F', 'PHM': 'F', 'PIV': 'X',
    'PLE': 'L', 'PM3': 'F', 'PMT': 'C', 'POM': 'P', 'PPN': 'F', 'PPU': 'A',
    'PPW': 'G', 'PQ1': 'N', 'PR3': 'C', 'PR5': 'A', 'PR9': 'P', 'PRN': 'A',
    'PRO': 'P', 'PRS': 'P', 'PSA': 'F', 'PSH': 'H', 'PST': 'T', 'PSU': 'U',
    'PSW': 'C', 'PTA': 'X', 'PTH': 'Y', 'PTM': 'Y', 'PTR': 'Y', 'PU': 'A',
    'PUY': 'N', 'PVH': 'H', 'PVL': 'X', 'PYA': 'A', 'PYO': 'U', 'PYX': 'C',
    'PYY': 'N', 'QMM': 'Q', 'QPA': 'C', 'QPH': 'F', 'QUO': 'G', 'R': 'A',
    'R1A': 'C', 'R4K': 'W', 'RE0': 'W', 'RE3': 'W', 'RIA': 'A', 'RMP': 'A',
    'RON': 'X', 'RT': 'T', 'RTP': 'N', 'S1H': 'S', 'S2C': 'C', 'S2D': 'A',
    'S2M': 'T', 'S2P': 'A', 'S4A': 'A', 'S4C': 'C', 'S4G': 'G', 'S4U': 'U',
    'S6G': 'G', 'SAC': 'S', 'SAH': 'C', 'SAR': 'G', 'SBL': 'S', 'SC': 'C',
    'SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SD2': 'X', 'SDG': 'G', 'SDP': 'S',
    'SEB': 'S', 'SEC': 'A', 'SEG': 'A', 'SEL': 'S', 'SEM': 'S', 'SEN': 'S',
    'SEP': 'S', 'SER': 'S', 'SET': 'S', 'SGB': 'S', 'SHC': 'C', 'SHP': 'G',
    'SHR': 'K', 'SIB': 'C', 'SLA': 'P', 'SLR': 'P', 'SLZ': 'K', 'SMC': 'C',
    'SME': 'M', 'SMF': 'F', 'SMP': 'A', 'SMT': 'T', 'SNC': 'C', 'SNN': 'N',
    'SOC': 'C', 'SOS': 'N', 'SOY': 'S', 'SPT': 'T', 'SRA': 'A', 'SSU': 'U',
    'STY': 'Y', 'SUB': 'X', 'SUN': 'S', 'SUR': 'U', 'SVA': 'S', 'SVV': 'S',
    'SVW': 'S', 'SVX': 'S', 'SVY': 'S', 'SVZ': 'X', 'SYS': 'C', 'T': 'T',
    'T11': 'F', 'T23': 'T', 'T2S': 'T', 'T2T': 'N', 'T31': 'U', 'T32': 'T',
    'T36': 'T', 'T37': 'T', 'T38': 'T', 'T39': 'T', 'T3P': 'T', 'T41': 'T',
    'T48': 'T', 'T49': 'T', 'T4S': 'T', 'T5O': 'U', 'T5S': 'T', 'T66': 'X',
    'T6A': 'A', 'TA3': 'T', 'TA4': 'X', 'TAF': 'T', 'TAL': 'N', 'TAV': 'D',
    'TBG': 'V', 'TBM': 'T', 'TC1': 'C', 'TCP': 'T', 'TCQ': 'Y', 'TCR': 'W',
    'TCY': 'A', 'TDD': 'L', 'TDY': 'T', 'TFE': 'T', 'TFO': 'A', 'TFQ': 'F',
    'TFT': 'T', 'TGP': 'G', 'TH6': 'T', 'THC': 'T', 'THO': 'X', 'THR': 'T',
    'THX': 'N', 'THZ': 'R', 'TIH': 'A', 'TLB': 'N', 'TLC': 'T', 'TLN': 'U',
    'TMB': 'T', 'TMD': 'T', 'TNB': 'C', 'TNR': 'S', 'TOX': 'W', 'TP1': 'T',
    'TPC': 'C', 'TPG': 'G', 'TPH': 'X', 'TPL': 'W', 'TPO': 'T', 'TPQ': 'Y',
    'TQI': 'W', 'TQQ': 'W', 'TRF': 'W', 'TRG': 'K', 'TRN': 'W', 'TRO': 'W',
    'TRP': 'W', 'TRQ': 'W', 'TRW': 'W', 'TRX': 'W', 'TS': 'N', 'TST': 'X',
    'TT': 'N', 'TTD': 'T', 'TTI': 'U', 'TTM': 'T', 'TTQ': 'W', 'TTS': 'Y',
    'TY1': 'Y', 'TY2': 'Y', 'TY3': 'Y', 'TY5': 'Y', 'TYB': 'Y', 'TYI': 'Y',
    'TYJ': 'Y', 'TYN': 'Y', 'TYO': 'Y', 'TYQ': 'Y', 'TYR': 'Y', 'TYS': 'Y',
    'TYT': 'Y', 'TYU': 'N', 'TYW': 'Y', 'TYX': 'X', 'TYY': 'Y', 'TZB': 'X',
    'TZO': 'X', 'U': 'U', 'U25': 'U', 'U2L': 'U', 'U2N': 'U', 'U2P': 'U',
    'U31': 'U', 'U33': 'U', 'U34': 'U', 'U36': 'U', 'U37': 'U', 'U8U': 'U',
    'UAR': 'U', 'UCL': 'U', 'UD5': 'U', 'UDP': 'N', 'UFP': 'N', 'UFR': 'U',
    'UFT': 'U', 'UMA': 'A', 'UMP': 'U', 'UMS': 'U', 'UN1': 'X', 'UN2': 'X',
    'UNK': 'X', 'UR3': 'U', 'URD': 'U', 'US1': 'U', 'US2': 'U', 'US3': 'T',
    'US5': 'U', 'USM': 'U', 'VAD': 'V', 'VAF': 'V', 'VAL': 'V', 'VB1': 'K',
    'VDL': 'X', 'VLL': 'X', 'VLM': 'X', 'VMS': 'X', 'VOL': 'X', 'X': 'G',
    'X2W': 'E', 'X4A': 'N', 'XAD': 'A', 'XAE': 'N', 'XAL': 'A', 'XAR': 'N',
    'XCL': 'C', 'XCN': 'C', 'XCP': 'X', 'XCR': 'C', 'XCS': 'N', 'XCT': 'C',
    'XCY': 'C', 'XGA': 'N', 'XGL': 'G', 'XGR': 'G', 'XGU': 'G', 'XPR': 'P',
    'XSN': 'N', 'XTH': 'T', 'XTL': 'T', 'XTR': 'T', 'XTS': 'G', 'XTY': 'N',
    'XUA': 'A', 'XUG': 'G', 'XX1': 'K', 'Y': 'A', 'YCM': 'C', 'YG': 'G',
    'YOF': 'Y', 'YRR': 'N', 'YYG': 'G', 'Z': 'C', 'Z01': 'A', 'ZAD': 'A',
    'ZAL': 'A', 'ZBC': 'C', 'ZBU': 'U', 'ZCL': 'F', 'ZCY': 'C', 'ZDU': 'U',
    'ZFB': 'X', 'ZGU': 'G', 'ZHP': 'N', 'ZTH': 'T', 'ZU0': 'T', 'ZZJ': 'A',
}
# common_typos_enable
# pyformat: enable


@functools.lru_cache(maxsize=64)
def letters_three_to_one(restype: str, *, default: str) -> str:
  """Returns single letter name if one exists otherwise returns default."""
  return CCD_NAME_TO_ONE_LETTER.get(restype, default)


ALA = sys.intern('ALA')
ARG = sys.intern('ARG')
ASN = sys.intern('ASN')
ASP = sys.intern('ASP')
CYS = sys.intern('CYS')
GLN = sys.intern('GLN')
GLU = sys.intern('GLU')
GLY = sys.intern('GLY')
HIS = sys.intern('HIS')
ILE = sys.intern('ILE')
LEU = sys.intern('LEU')
LYS = sys.intern('LYS')
MET = sys.intern('MET')
PHE = sys.intern('PHE')
PRO = sys.intern('PRO')
SER = sys.intern('SER')
THR = sys.intern('THR')
TRP = sys.intern('TRP')
TYR = sys.intern('TYR')
VAL = sys.intern('VAL')
UNK = sys.intern('UNK')
GAP = sys.intern('-')

# Unknown ligand.
UNL = sys.intern('UNL')

# Non-standard version of MET (with Se instead of S), but often appears in PDB.
MSE = sys.intern('MSE')

# 20 standard protein amino acids (no unknown).
PROTEIN_TYPES: tuple[str, ...] = (
    ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
)  # pyformat: disable

# 20 standard protein amino acids plus the unknown (UNK) amino acid.
PROTEIN_TYPES_WITH_UNKNOWN: tuple[str, ...] = PROTEIN_TYPES + (UNK,)

# This is the standard residue order when coding AA type as a number.
# Reproduce it by taking 3-letter AA codes and sorting them alphabetically.
# For legacy reasons this only refers to protein residues.

PROTEIN_TYPES_ONE_LETTER: tuple[str, ...] = (
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V',
)  # pyformat: disable

PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN: tuple[str, ...] = (
    PROTEIN_TYPES_ONE_LETTER + ('X',)
)
PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_AND_GAP: tuple[str, ...] = (
    PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN + (GAP,)
)

PROTEIN_TYPES_ONE_LETTER_TO_INT: Mapping[str, int] = {
    r: i for i, r in enumerate(PROTEIN_TYPES_ONE_LETTER)
}
PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_TO_INT: Mapping[str, int] = {
    r: i for i, r in enumerate(PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN)
}

PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_AND_GAP_TO_INT: Mapping[str, int] = {
    r: i for i, r in enumerate(PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_AND_GAP)
}


PROTEIN_COMMON_ONE_TO_THREE: Mapping[str, str] = {
    'A': ALA,
    'R': ARG,
    'N': ASN,
    'D': ASP,
    'C': CYS,
    'Q': GLN,
    'E': GLU,
    'G': GLY,
    'H': HIS,
    'I': ILE,
    'L': LEU,
    'K': LYS,
    'M': MET,
    'F': PHE,
    'P': PRO,
    'S': SER,
    'T': THR,
    'W': TRP,
    'Y': TYR,
    'V': VAL,
}

PROTEIN_COMMON_THREE_TO_ONE: Mapping[str, str] = {
    v: k for k, v in PROTEIN_COMMON_ONE_TO_THREE.items()
}

A = sys.intern('A')
G = sys.intern('G')
C = sys.intern('C')
U = sys.intern('U')
T = sys.intern('T')

DA = sys.intern('DA')
DG = sys.intern('DG')
DC = sys.intern('DC')
DT = sys.intern('DT')

UNK_NUCLEIC_ONE_LETTER = sys.intern('N')  # Unknown nucleic acid single letter.
UNK_RNA = sys.intern('N')  # Unknown RNA.
UNK_DNA = sys.intern('DN')  # Unknown DNA residue (differs from N).

RNA_TYPES: tuple[str, ...] = (A, G, C, U)
DNA_TYPES: tuple[str, ...] = (DA, DG, DC, DT)

NUCLEIC_TYPES: tuple[str, ...] = RNA_TYPES + DNA_TYPES
# Without UNK DNA.
NUCLEIC_TYPES_WITH_UNKNOWN: tuple[str, ...] = NUCLEIC_TYPES + (
    UNK_NUCLEIC_ONE_LETTER,
)
NUCLEIC_TYPES_WITH_2_UNKS: tuple[str, ...] = NUCLEIC_TYPES + (
    UNK_RNA,
    UNK_DNA,
)

RNA_TYPES_ONE_LETTER_WITH_UNKNOWN: tuple[str, ...] = RNA_TYPES + (UNK_RNA,)
RNA_TYPES_ONE_LETTER_WITH_UNKNOWN_TO_INT: Mapping[str, int] = {
    r: i for i, r in enumerate(RNA_TYPES_ONE_LETTER_WITH_UNKNOWN)
}

DNA_TYPES_WITH_UNKNOWN: tuple[str, ...] = DNA_TYPES + (UNK_DNA,)
DNA_TYPES_ONE_LETTER: tuple[str, ...] = (A, G, C, T)
DNA_TYPES_ONE_LETTER_WITH_UNKNOWN: tuple[str, ...] = DNA_TYPES_ONE_LETTER + (
    UNK_NUCLEIC_ONE_LETTER,
)
DNA_TYPES_ONE_LETTER_WITH_UNKNOWN_TO_INT: Mapping[str, int] = {
    r: i for i, r in enumerate(DNA_TYPES_ONE_LETTER_WITH_UNKNOWN)
}
DNA_COMMON_ONE_TO_TWO: Mapping[str, str] = {
    'A': 'DA',
    'G': 'DG',
    'C': 'DC',
    'T': 'DT',
}

STANDARD_POLYMER_TYPES: tuple[str, ...] = PROTEIN_TYPES + NUCLEIC_TYPES
POLYMER_TYPES: tuple[str, ...] = PROTEIN_TYPES_WITH_UNKNOWN + NUCLEIC_TYPES
POLYMER_TYPES_WITH_UNKNOWN: tuple[str, ...] = (
    PROTEIN_TYPES_WITH_UNKNOWN + NUCLEIC_TYPES_WITH_UNKNOWN
)
POLYMER_TYPES_WITH_GAP: tuple[str, ...] = PROTEIN_TYPES + (GAP,) + NUCLEIC_TYPES
POLYMER_TYPES_WITH_UNKNOWN_AND_GAP: tuple[str, ...] = (
    PROTEIN_TYPES_WITH_UNKNOWN + (GAP,) + NUCLEIC_TYPES_WITH_UNKNOWN
)
POLYMER_TYPES_WITH_ALL_UNKS_AND_GAP: tuple[str, ...] = (
    PROTEIN_TYPES_WITH_UNKNOWN + (GAP,) + NUCLEIC_TYPES_WITH_2_UNKS
)

POLYMER_TYPES_ORDER = {restype: i for i, restype in enumerate(POLYMER_TYPES)}

POLYMER_TYPES_ORDER_WITH_UNKNOWN = {
    restype: i for i, restype in enumerate(POLYMER_TYPES_WITH_UNKNOWN)
}

POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP = {
    restype: i for i, restype in enumerate(POLYMER_TYPES_WITH_UNKNOWN_AND_GAP)
}

POLYMER_TYPES_ORDER_WITH_ALL_UNKS_AND_GAP = {
    restype: i for i, restype in enumerate(POLYMER_TYPES_WITH_ALL_UNKS_AND_GAP)
}

POLYMER_TYPES_NUM = len(POLYMER_TYPES)  # := 29.
POLYMER_TYPES_NUM_WITH_UNKNOWN = len(POLYMER_TYPES_WITH_UNKNOWN)  # := 30.
POLYMER_TYPES_NUM_WITH_GAP = len(POLYMER_TYPES_WITH_GAP)  # := 29.
POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP = len(
    POLYMER_TYPES_WITH_UNKNOWN_AND_GAP
)  # := 31.
POLYMER_TYPES_NUM_ORDER_WITH_ALL_UNKS_AND_GAP = len(
    POLYMER_TYPES_WITH_ALL_UNKS_AND_GAP
)  # := 32.

WATER_TYPES: tuple[str, ...] = ('HOH', 'DOD')

UNKNOWN_TYPES: tuple[str, ...] = (UNK, UNK_RNA, UNK_DNA, UNL)
