# Known Issues

## Numerical performance for CUDA Capability 7.x GPUs

All CUDA Capability 7.x GPUs (e.g. V100) produce obviously bad output, with lots
of clashing residues (the clashes cause a ranking score of -99 or lower), unless
the environment variable `XLA_FLAGS` is set to include
`--xla_disable_hlo_passes=custom-kernel-fusion-rewriter`.

## Incorrect handling of two-letter atoms in SMILES ligands

Between commits https://github.com/google-deepmind/alphafold3/commit/f8df1c7 and
https://github.com/google-deepmind/alphafold3/commit/4e4023c, AlphaFold 3
handled incorrectly any two-letter atoms (e.g. Cl, Br) in ligands defined using
SMILES strings.

**Current status:**  
This issue was resolved starting from commit 4e4023c (which standardized atom name assignment using a shared RDKit utility). Two-letter elements like Cl and Br in SMILES are now parsed and handled correctly on the main branch.

To demonstrate and validate proper input for common drug-like molecules (halogens appear frequently in pharmaceuticals), a new example has been added:

- File: `test_inputs/ligand_halogen_cl.json`  
- SMILES: `CC(=O)Nc1ccc(Cl)cc1` (paracetamol analog with chlorine)

This serves as a reference test case for users working with halogen-containing ligands.

## MSA discrepancy between AlphaFold 3 and AlphaFold Server

### The root cause of the problem

The released AlphaFold 3 and AlphaFold Server use the same model weights and
equivalent featurisation and model code. However, the way they run genetic
search is slightly different. The released AlphaFold 3 searches each database in
one go, while AlphaFold Server has a sharded version of each database (split
into multiple smaller FASTA files) and searches all of the shards in parallel.
The results of these parallel searches are then merged together at the end.

The discrepancy is caused by a different (deeper) MSA on AlphaFold Server in
some cases. We discovered that the issue is caused by running sharded Jackhmmer
in AlphaFold Server without the `--domZ` flag (has to be set together with the
`--Z` flag and set to the same value) which means that effectively the AlphaFold
Server is running with roughly 100× more permissive `--domE` filter. This means
more sequences are sometimes included in the MSA.

We are keeping behaviour unchanged in both the released AlphaFold 3 and in the
AlphaFold Server, however, we are giving users with local installs an option to
replicate AlphaFold Server behaviour locally. In our large scale tests the
difference did not matter, it is only very specific inputs that get better
accuracy with the deeper MSA.

See https://github.com/google-deepmind/alphafold3/issues/492 for an example
input where a protein-DNA complex gets significantly higher ipTM and pTM with
AlphaFold Server compared to a local run.

### Replicating AlphaFold Server behaviour locally

If you want to replicate AlphaFold Server behaviour (i.e. better folding
accuracy in some cases), you can increase the value of the Jackhmmer/Nhmmer
`--domE` flag by 100× compared to its default value.

Alternatively, you can run the sharded MSA search while not setting the `--domZ`
value – you would have to modify the code to do it. We added support for
searching against sharded databases in AlphaFold 3 in
https://github.com/google-deepmind/alphafold3/commit/805adc3863841d83d631ccd18136ad58ce3ecb34
and the way to run AlphaFold 3 with sharded databases is documented in
https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#sharded-genetic-databases.
It can provide 10–30× speedup (potentially even more, depending on hardware) of
the genetic search.

In general, we recommend experimenting with MSA if you are seeing a prediction
with low predicted confidence. Typically adding more *relevant* sequences in the
MSA will increase AlphaFold prediction accuracy and model confidence scores.
