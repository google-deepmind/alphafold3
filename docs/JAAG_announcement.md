# JAAG — JSON input file Assembler for AlphaFold 3 (with Glycan integration)

Note: JAAG is fully compatible with standalone AlphaFold 3 but is not
compatible with the AlphaFold 3 public server.

We’re excited to share JAAG (JSON input file Assembler for AlphaFold 3) — a
lightweight web GUI that helps researchers build AlphaFold 3 input `.json`
files with correct glycan syntax (including bondedAtomPairs + CCD formatting)
so glycans get improved stereochemistry during AF3 inference.

JAAG automates composing AF3 input JSONs for glycoproteins and glycan–protein
complexes so users don't have to hand-edit the often error-prone glycan
sections. The web app generates a complete AF3 input file containing the
expected configuration fields and glycan-specific blocks.

Resources
- Web app: https://biofgreat.org/JAAG
- Source code (open-source): https://github.com/chinchc/JAAG
- Preprint: https://www.biorxiv.org/content/10.1101/2025.10.08.681238v1
- Reference on interpreting AF3 glyco-related models: https://doi.org/10.1093/glycob/cwaf048

Privacy
No user data or JSON content is stored or processed on the JAAG server.

Contact
Warm regards,

Chin Huang | Moremen Lab
Ph.D. Candidate & Graduate Research Assistant
Department of Biochemistry and Molecular Biology
Complex Carbohydrate Research Center
University of Georgia

---

If you would like, I can:
- Prepare a polished GitHub issue body for maintainers and open the issue.
- Create a pull request that adds this announcement to `docs/` (I can create a
  branch, add this file, run basic checks, and open the PR).
- Produce three polished variants of this message (issue, PR description, and
  mailing-list/post format) so you can copy them as needed.
