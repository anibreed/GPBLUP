# GPBLUP — synthetic demonstration dataset

This folder lets anyone **run GPBLUP end to end and reproduce a representative
single-step genomic evaluation without any real animal data**. The individual-level
pig genotype and pedigree records used to validate GPBLUP in the paper are owned by a
commercial breeding company and cannot be redistributed; the data here are **fully
synthetic** and carry no real-animal information.

## What the demo does

A small pig-like nucleus is simulated: 60 founders plus four generations of 160
animals (700 total), with genotypes **gene-dropped through the pedigree** so that
relatives share genotypes and the genomic relationship matrix reflects the pedigree.
Two correlated traits (ADG, BF) are generated from marker effects plus fixed
(herd-year, sex), covariate (age), litter (IID) and residual terms. The last three
generations (480 animals) are genotyped; the first four generations (640 animals) are
phenotyped. GPBLUP then runs a two-trait single-step (H) evaluation with the APY
engine, and `eval_demo.py` correlates the resulting GEBV with the hidden true breeding
values. Everything is deterministic (fixed seed 2026).

- 700 animals (480 genotyped, 640 phenotyped) × 800 SNP, 2 traits.

## Files

| File | Role |
|------|------|
| `ped.txt` | pedigree: `animal sire dam` (0 = unknown parent) |
| `geno.txt` | genotypes: `animal dosage-string` (one 0/1/2 char per SNP) |
| `phen.txt` | phenotypes (header): `animal dummy herd_year sex ADG BF age litter` |
| `model_demo.par` | the parameter file that drives the run |
| `demo_truth.txt` | hidden true breeding values (for the accuracy check) |
| `make_demo.py` | regenerates all of the above (needs NumPy) |
| `eval_demo.py` | scores GEBV against the true breeding values |
| `run_demo.sh` | one-command runner |

## Run

```bash
chmod +x ../bin/*
./run_demo.sh
```

## Expected result

GPBLUP's PCG solver converges, writes GEBV to `phen.txt.ebv`, and the fixed-effect
solutions recover the simulated values (e.g., the age covariate ≈ −0.20 / 0.05 and the
sex contrast ≈ −2 / +1). `eval_demo.py` reports the GEBV-versus-true-breeding-value
correlations, which are clearly positive for both traits (of order 0.6–0.7 on this
small synthetic set), confirming the single-step evaluation works end to end.

Generated files (`phen.txt.ebv`, `*_ren`, `*.ren.idx`, `*_inbcache`) are recreated on
every run and are not tracked.
