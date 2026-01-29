# Community Tools

## JAAG: a JSON input file Assembler for AlphaFold 3 (with Glycan Integration)

JAAG is a lightweight, web-based GUI tool that helps generate AlphaFold 3 input
JSON files with integrated glycan support. It automates the creation of correct
glycan syntax (including `bondedAtomPairs` + CCD), reducing manual errors when
preparing glycoprotein or glycan–protein complexes.

*   Web app: https://biofgreat.org/JAAG
*   Source code: https://github.com/chinchc/JAAG
*   Paper: https://doi.org/10.1093/glycob/cwaf083

Note: JAAG is compatible with standalone AlphaFold 3, but not with the AlphaFold
3 server.

## Modeling glycans with AlphaFold 3: capabilities, caveats, and limitations

Paper on modeling glycans (and other ligands) with AF3 that modeled and assessed
major glycan classes and provides:

*   Step-by-step tutorial for building ligand inputs (applicable beyond glycans)
*   Ready-to-run scripts for each glycan class
*   Comprehensive CCD table for all SNFG monosaccharides
*   Discussion of caveats and limitations of AF3
*   Full AF3 inputs/outputs archived on ModelArchive for reproducibility

Useful resource if your AF3 ligand models appear stereochemically off.

*   Paper: https://doi.org/10.1093/glycob/cwaf048
*   ModelArchive: https://doi.org/10.5452/ma-af3glycan
