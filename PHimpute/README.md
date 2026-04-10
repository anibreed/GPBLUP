# PHimpute

Missing SNP imputation from valid SNP datasets.

## Purpose

- Impute missing SNP genotypes after population QC.
- Preserve numeric genotype encoding and missing value conventions.

## Usage

```bash
phimpute <parameter_file>
```

See USER_GUIDE.md in the repository root for the pipeline context.

## Source Layout
- PHimpute main and PHimpute-only modules live under PHimpute/ and PHimpute/src.
- Shared core modules remain under src/ (used by multiple programs).

## Module Boundaries
- PHimpute uses shared modules only when the functionality is generic and reused elsewhere.
- PHimpute-only modules stay under PHimpute/src and do not import other programs.

## Shared Module Dependencies
PHimpute currently relies on these shared modules from src/:

- M_Kinds: numeric kinds shared across programs.
- M_Variables: shared FileInfo structure and species configuration utilities.
- M_ReadFile: file open and record reading helpers.
- M_StrEdit: string parsing utilities (delimiters, tokenizing, case conversion).
- M_IO_Ped: pedigree IO helpers used by deterministic pedigree imputation.
- M_ped: pedigree record types and hash utilities.

If a dependency becomes PHimpute-specific, move the relevant code into PHimpute/src
and update imports to keep shared modules minimal.
