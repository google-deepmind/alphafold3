# AlphaFold 3 Output

## Output Directory Structure

For every input job, AlphaFold 3 writes all its outputs in a directory called by
the sanitized version of the job name. E.g. for job name "My first fold (test)",
AlphaFold 3 will write its outputs in a directory called `my_first_fold_test`.
If such directory already exists, AlphaFold 3 will append a timestamp to the
directory name to avoid overwriting existing data.

The following structure is used within the output directory:

*   Sub-directories with results for each sample and seed. There will be
    *num\_seeds* \* *num\_samples* such sub-directories. The naming pattern is
    `seed-<seed value>_sample-<sample number>`. Each of these directories
    contains a confidence JSON, summary confidence JSON, and the mmCIF with the
    predicted structure.
*   Embeddings for each seed: `seed-<seed value>_embeddings/embeddings.npz`.
    Only saved if AlphaFold 3 is run with `--save_embeddings=true`.
*   Top-ranking prediction mmCIF: `<job_name>_model.cif`. This file contains the
    predicted coordinates and should be compatible with most structural biology
    tools. We do not provide the output in the PDB format, the CIF file can be
    easily converted into one if needed.
*   Top-ranking prediction confidence JSON: `<job_name>_confidences.json`.
*   Top-ranking prediction summary confidence JSON:
    `<job_name>_summary_confidences.json`.
*   Job input JSON file with the MSA and template data added by the data
    pipeline: `<job_name>_data.json`.
*   Ranking scores for all predictions: `ranking_scores.csv`. The prediction
    with highest ranking is the one included in the root directory.
*   Output terms of use: `TERMS_OF_USE.md`.

Below is an example AlphaFold 3 output directory listing for a job called "Hello
Fold", that has been ran with 1 seed and 5 samples:

```txt
hello_fold/
├── seed-1234_embeddings                       # Only if --save_embeddings=true.
│   └── hello_fold_seed-1234_embeddings.npz    # Only if --save_embeddings=true.
├── seed-1234_sample-0/
│   ├── hello_fold_seed-1234_sample-0_confidences.json
│   ├── hello_fold_seed-1234_sample-0_model.cif
│   └── hello_fold_seed-1234_sample-0_summary_confidences.json
├── seed-1234_sample-1/
│   ├── hello_fold_seed-1234_sample-1_confidences.json
│   ├── hello_fold_seed-1234_sample-1_model.cif
│   └── hello_fold_seed-1234_sample-1_summary_confidences.json
├── seed-1234_sample-2/
│   ├── hello_fold_seed-1234_sample-2_confidences.json
│   ├── hello_fold_seed-1234_sample-2_model.cif
│   └── hello_fold_seed-1234_sample-2_summary_confidences.json
├── seed-1234_sample-3/
│   ├── hello_fold_seed-1234_sample-3_confidences.json
│   ├── hello_fold_seed-1234_sample-3_model.cif
│   └── hello_fold_seed-1234_sample-3_summary_confidences.json
├── seed-1234_sample-4/
│   ├── hello_fold_seed-1234_sample-4_confidences.json
│   ├── hello_fold_seed-1234_sample-4_model.cif
│   └── hello_fold_seed-1234_sample-4_summary_confidences.json
├── TERMS_OF_USE.md
├── hello_fold_confidences.json
├── hello_fold_data.json
├── hello_fold_model.cif
├── hello_fold_ranking_scores.csv
└── hello_fold_summary_confidences.json
```

## Confidence Metrics

Similar to AlphaFold 2 and AlphaFold-Multimer, AlphaFold 3 outputs include
confidence metrics. The main metrics are:

*   **pLDDT:** a per-atom confidence estimate on a 0-100 scale where a higher
    value indicates higher confidence. pLDDT aims to predict a modified LDDT
    score that only considers distances to polymers. For proteins this is
    similar to the
    [lDDT-Cα metric](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3799472/) but
    with more granularity as it can vary per atom not just per residue. For
    ligand atoms, the modified LDDT considers the errors only between the ligand
    atom and polymers, not other ligand atoms. For DNA/RNA a wider radius of 30
    Å is used for the modified LDDT instead of 15 Å.
*   **PAE (predicted aligned error)**: an estimate of the error in the relative
    position and orientation between two tokens in the predicted structure.
    Higher values indicate higher predicted error and therefore lower
    confidence. For proteins and nucleic acids, PAE score is essentially the
    same as AlphaFold 2, where the error is measured relative to frames
    constructed from the protein backbone. For small molecules and
    post-translational modifications, a frame is constructed for each atom from
    its closest neighbors from a reference conformer.
*   **pTM and ipTM scores**: the predicted template modeling (pTM) score and the
    interface predicted template modeling (ipTM) score are both derived from a
    measure called the template modeling (TM) score. This measures the accuracy
    of the entire structure
    ([Zhang and Skolnick, 2004](https://doi.org/10.1002/prot.20264);
    [Xu and Zhang, 2010](https://doi.org/10.1093/bioinformatics/btq066)). A pTM
    score above 0.5 means the overall predicted fold for the complex might be
    similar to the true structure. ipTM measures the accuracy of the predicted
    relative positions of the subunits within the complex. Values higher than
    0.8 represent confident high-quality predictions, while values below 0.6
    suggest a failed prediction. ipTM values between 0.6 and 0.8 are a gray zone
    where predictions could be correct or incorrect. The TM score is very strict
    for small structures or short chains, so pTM assigns values less than 0.05
    when fewer than 20 tokens are involved; for these cases PAE or pLDDT may be
    more indicative of prediction quality.

For detailed description of these confidence metrics see the
[AlphaFold 3 paper](https://www.nature.com/articles/s41586-024-07487-w). For
protein components, the
[AlphaFold: A Practical guide](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/)
course for structures provides additional tutorials on the confidence metrics.

If you are interested in a specific entity or interaction, then there are
confidences available in the outputs which are specific to each chain or
chain-pair, as opposed to the full complex. See below for more details on all
the confidence metrics that are returned.

## Multi-Seed and Multi-Sample Results

By default, the model samples five predictions per seed. The top-ranked
prediction across all samples and seeds is available at the top-level of the
output directory. All samples along with their associated confidences are
available in subdirectories of the output directory.

For ranking of the full complex use the `ranking_score` (higher is better). This
score uses overall structure confidences (pTM and ipTM), but also includes terms
that penalize clashes and encourage disordered regions not to have spurious
helices – these extra terms mean the score should only be used to rank
structures.

If you are interested in a specific entity or interaction, you may want to rank
by a metric specific to that chain or chain-pair, as opposed to the full
complex. In that case, use the per chain or per chain-pair confidence metrics
described below for ranking.

## Metrics in Confidences JSON

For each predicted sample we provide two JSON files. One contains summary
metrics – summaries for either the whole structure, per chain or per chain-pair
– and the other contains full 1D or 2D arrays.

Summary outputs:

*   `ptm`: A scalar in the range 0-1 indicating the predicted TM-score for the
    full structure.
*   `iptm`: A scalar in the range 0-1 indicating predicted interface TM-score
    (confidence in the predicted interfaces) for all interfaces in the
    structure.
*   `fraction_disordered`: A scalar in the range 0-1 that indicates what
    fraction of the prediction structure is disordered, as measured by
    accessible surface area, see our
    [paper](https://www.nature.com/articles/s41586-024-07487-w) for details.
*   `has_clash`: A boolean indicating if the structure has a significant number
    of clashing atoms (more than 50% of a chain, or a chain with more than 100
    clashing atoms).
*   `ranking_score`: A scalar in the range \[-100, 1.5\] that can be used for
    ranking predictions, it incorporates `ptm`, `iptm`, `fraction_disordered`
    and `has_clash` into a single number with the following equation: 0.8 × ipTM
    \+ 0.2 × pTM \+ 0.5 × disorder − 100 × has_clash.
*   `chain_pair_pae_min`: A \[num_chains, num_chains\] array. Element (i, j) of
    the array contains the lowest PAE value across rows restricted to chain i
    and columns restricted to chain j. This has been found to correlate with
    whether two chains interact or not, and in some cases can be used to
    distinguish binders from non-binders.
*   `chain_pair_iptm`: A \[num_chains, num_chains\] array. Off-diagonal element
    (i, j) of the array contains the ipTM restricted to tokens from chains i and
    j. Diagonal element (i, i) contains the pTM restricted to chain i. Can be
    used for ranking a specific interface between two chains, when you know that
    they interact, e.g. for antibody-antigen interactions
*   `chain_ptm`: A \[num_chains\] array. Element i contains the pTM restricted
    to chain i. Can be used for ranking individual chains when the structure of
    that chain is most of interest, rather than the cross-chain interactions it
    is involved with.
*   `chain_iptm:` A \[num_chains\] array that gives the average confidence
    (interface pTM) in the interface between each chain and all other chains.
    Can be used for ranking a specific chain, when you care about where the
    chain binds to the rest of the complex and you do not know which other
    chains you expect it to interact with. This is often the case with ligands.

Full array outputs:

*   `pae`: A \[num\_tokens, num\_tokens\] array. Element (i, j) indicates the
    predicted error in the position of token j, when the prediction is aligned
    to the ground truth using the frame of token i.
*   `atom_plddts`: A \[num_atoms\] array, element i indicates the predicted
    local distance difference test (pLDDT) for atom i in the prediction.
*   `contact_probs`: A \[num_tokens, num_tokens\] array. Element (i, j)
    indicates the predicted probability that token i and token j are in contact
    (8 Å between the representative atom for each token), see
    [paper](https://www.nature.com/articles/s41586-024-07487-w) for details.
*   `token_chain_ids`: A \[num_tokens\] array indicating the chain ids
    corresponding to each token in the prediction.
*   `atom_chain_ids`: A \[num_atoms\] array indicating the chain ids
    corresponding to each atom in the prediction.

## Embeddings

AlphaFold 3 can be run with `--save_embeddings=true` to save the embeddings for
each seed. The file is in the
[compressed Numpy `.npz` format](https://numpy.org/doc/stable/reference/generated/numpy.savez_compressed.html)
and can be loaded using `numpy.load` as a dictionary-like object with two
arrays:

*   `single_embeddings`: A \`[num\_tokens, 384\] array containing the embeddings
    for each token.
*   `pair_embeddings`: A \[num\_tokens, num\_tokens, 128\] array containing the
    pairwise embeddings between all tokens.

You can use for instance the following Python code to load the embeddings:

```py
import numpy as np

with open('embeddings.npz', 'rb') as f:
  embeddings = np.load(f)
  single_embeddings = embeddings['single_embeddings']
  pair_embeddings = embeddings['pair_embeddings']
```

## Chirality checks

In the AlphaFold 3 paper Posebusters results, a penalty was applied to the
ranking score if the ligand of interest contained chiral errors. By running
multiple seeds and using this chiral aware ranking, chiral error rates were
greatly reduced.

We provide the method `compare_chirality` in
[`model/scoring/chirality.py`](https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/scoring/chirality.py)
to replicate these chiral checks. Chirality is checked against CCD structures if
available, otherwise users can supply custom RDKit Mol objects for comparison.
