# Algorithm Background

This document summarizes the algorithms and data flow across the GPBLUP toolchain.
It focuses on the high-level logic users need to understand the pipeline.

## Pipeline Overview

Input data (Illumina FinalReport, PED, MAP)
  -> ReadFR  (numeric conversion + initial QC)
  -> popQC   (population-based SNP/animal QC)
  -> PHimpute (missing SNP imputation)
  -> relgeno (genomic relationship matrix G and G inverse)
  -> relped  (pedigree relationship matrix A and A inverse)

## ReadFR (FinalReport to numeric GENO)

Purpose:
- Convert Illumina FinalReport records to numeric SNP genotypes.
- Apply per-record QC thresholds and output a GENO matrix.

Key steps:
- Parse PED, MAP, and FinalReport columns using a parameter file.
- Map sample identifiers (ANIMAL_ID or ANIMAL_ARN) to pedigree records.
- Filter SNP records by quality thresholds (GC score, R intensity, GT score,
  cluster separation, call rate).
- Encode genotype as 0/1/2 and emit missing as a dedicated code (e.g., 9).

Outputs:
- GENO file with numeric SNP genotypes and pedigree context.

## popQC (Population QC)

Purpose:
- Evaluate and filter SNPs and animals at the population level.

SNP-level QC (typical):
- Call rate: remove SNPs with high missingness.
- MAF: remove rare variants below a threshold.
- HWE: remove SNPs failing Hardy-Weinberg equilibrium tests.
- Monomorphic detection: remove non-informative SNPs.

Animal-level QC (typical):
- Sex concordance: compare X-chromosome heterozygosity to PED sex.
- Mendelian inconsistency: remove or flag parent-offspring conflicts.

Outputs:
- Filtered genotype data and QC reports.

## PHimpute (Missing SNP imputation)

Purpose:
- Impute missing SNP genotypes in the valid SNP set.

Typical approach:
- Use pedigree constraints when available.
- Use local SNP patterns to infer missing values.
- Preserve original genotype encoding and missing code conventions.

Outputs:
- Imputed genotype dataset used for downstream relationship matrices.

## relgeno (Genomic relationship matrix)

Purpose:
- Compute G and (optionally) G inverse from genotype data.

Common formulation (VanRaden-style):
- Let M be the genotype matrix (0,1,2) and p be allele frequencies.
- Centered matrix: Z = M - 2p
- Scaling: c = 2 * sum_j p_j (1 - p_j)
- G = (Z * Z^T) / c

Notes:
- A diagonal blend or weight can be applied for numerical stability.
- G inverse is computed if requested by parameters.

## relped (Pedigree relationship matrix)

Purpose:
- Compute A and (optionally) A inverse from pedigree records.

Typical method:
- Renumber pedigree records to satisfy parent-before-offspring order.
- Use Henderson's rules to build A inverse directly or via L D L^T factors.
- Output A and A inverse or selected formats for downstream use.

## Data Flow Notes

- Consistent sample identifiers are essential across PED, MAP, and GENO.
- Missing genotype code must be kept consistent across programs.
- Parameter files define all thresholds and column mappings.

## Related Documents

- Build and install: BUILD.md, INSTALL.md
- Usage: USER_GUIDE.md
- Program-specific details are under each program folder.
