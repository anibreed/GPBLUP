# GPBLUP User Guide

This guide summarizes how to run each program. For full details, use the
program-specific documentation listed below.

## End-to-End Workflow

1. ReadFR: convert Illumina FinalReport to numeric GENO with initial QC
2. popQC: population-level SNP and animal QC
3. phimpute: impute missing SNPs in the valid SNP set
4. relgeno: compute G and G inverse from SNP data
5. relped: compute A and A inverse from pedigree data

## ReadFR

Purpose: parse Illumina FinalReport files and generate GENO output.

```bash
ReadFR <parameter_file>
```

Docs: [ReadFR/README.md](ReadFR/README.md), [ReadFR/READFR_USER_MANUAL.md](ReadFR/READFR_USER_MANUAL.md)

## popQC

Purpose: SNP and animal QC pipeline.

```bash
popQC <parameter_file>
```

Docs: [popQC/README.md](popQC/README.md), [popQC/PIPELINE_GUIDE.md](popQC/PIPELINE_GUIDE.md)

## relped

Purpose: pedigree renumbering and A/Ainv matrix outputs.

```bash
relped <parameter_file>
```

Docs: [relped/README.md](relped/README.md)

## relgeno

Purpose: genomic relationship matrix processing.

```bash
relgeno <parameter_file>
```

Docs: [relgeno/RELGENO_PARAMS.md](relgeno/RELGENO_PARAMS.md)

## phimpute

Purpose: missing SNP imputation with pedigree- and population-aware steps.

```bash
phimpute <parameter_file>
```

Docs: [PHimpute/README.md](PHimpute/README.md), [PHimpute/PARAMETER_GUIDE.md](PHimpute/PARAMETER_GUIDE.md)

### PHimpute Algorithm (Detailed)

PHimpute executes a staged pipeline. Each stage can be enabled/disabled in the
parameter file.

1) Input loading and normalization
- MAP is read first; if ARRAY_ALL/ARRAY_CHR are missing or invalid, PHimpute
	rebuilds them by CHR/POS sorting.
- Candidate GENO is loaded using field mapping (ANIMAL_ID, GENO).
- Optional reference GENO is loaded; otherwise candidate GENO is used as the
	reference panel.

2) Optional QC filtering
- Sample filter: remove animals with missing rate > QC_SAMPLE_MISS_MAX.
- SNP filter: remove SNPs with missing rate > QC_SNP_MISS_MAX or MAF < QC_MIN_MAF.
- Optional ambiguous SNP filter (A/T or C/G) when allele info is available.
- Removed samples and SNPs are logged with the QC_LOG_PREFIX.

3) Pedigree-deterministic imputation
- For each child with known sire/dam, missing genotypes are filled using Mendel
	deterministic rules when parents are informative.
- This step uses the pedigree hash table (PEDFILE) and operates per chromosome.

4) Family phasing (trio-based)
- Build a haplotype panel from the reference GENO using sliding windows.
- For each trio (child/sire/dam), missing child genotypes are imputed by Mendel
	phase rules and tallied by chromosome and SNP.
- If a reference panel is provided, PHimpute compares fill rates with and
	without the reference panel before running the final pass.

5) Haplotype matching and imputation
- A haplotype panel is built using window and shift sizes (HAP_TUNE can sweep
	candidate values).
- For each window, if the candidate has enough non-missing SNPs, PHimpute
	hashes the window and finds matching reference haplotypes.
- Missing SNPs in the window are filled by weighted majority vote from the
	top-k matching reference haplotypes.

6) Population LD imputation (optional)
- For each missing genotype, PHimpute aggregates genotype counts using LD
	windows to the left and right of the target SNP.
- The genotype with the highest aggregated count is chosen.
- Window size is controlled by LD_WINDOW and LD_WINDOW_MAX.

7) Output
- The final imputed GENO file is written to OUTPUT: GENO.

## Additional References

- Algorithm background: [docs/ALGORITHM_BACKGROUND.md](docs/ALGORITHM_BACKGROUND.md)
- Model proof (sketches): [docs/MODEL_PROOF.md](docs/MODEL_PROOF.md)