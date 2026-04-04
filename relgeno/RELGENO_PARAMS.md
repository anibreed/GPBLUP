# relgeno parameter schema

This document describes the relgeno parameter file schema used by the relgeno program.

## Format

- One key per line, tokens separated by spaces.
- Keys are case-insensitive and should end with a colon.
- Lines can include comments starting with `#`.
- Paths can be quoted with single or double quotes.
- Relative paths are resolved relative to the parameter file location.

## Required

- `GENOFILE:` Path to the genotype input file.

## Optional outputs

- `OUT_G:` Output path for G matrix (default: `GENOFILE.G`).
- `OUT_GINV:` Output path for G inverse (default: `GENOFILE.Gi`).
- `OUT_GDIAG:` Output path for diagonal summary (default: `GENOFILE.Gdiag`).

## Optional parameters

- `MISSING:` Missing genotype code (default: `5`).
- `WEIGHT:` Weight factor for G (default: `0.95`).
- `DIAG_WEIGHT:` Diagonal adjustment (default: `1 - WEIGHT`).
- `COMPUTE_G:` `TRUE` or `FALSE` (default: `TRUE`).
- `COMPUTE_GINV:` `TRUE` or `FALSE` (default: `TRUE`).

## Aliases

- `OUT_GFILE:` or `G_OUT:` for `OUT_G:`
- `OUT_GI:` or `GINV_OUT:` for `OUT_GINV:`
- `OUT_G_DIAG:` for `OUT_GDIAG:`
- `MISSING_CODE:` for `MISSING:`
- `G_WEIGHT:` for `WEIGHT:`
- `DIAGONAL_WEIGHT:` for `DIAG_WEIGHT:`

## Example (Template)

```
# relgeno parameter file template (placeholders only)
GENOFILE: /path/to/genotypes.txt
OUT_G: /path/to/output/geno.G
OUT_GINV: /path/to/output/geno.Gi
OUT_GDIAG: /path/to/output/geno.Gdiag
MISSING: 5
WEIGHT: 0.95
DIAG_WEIGHT: 0.05
COMPUTE_G: TRUE
COMPUTE_GINV: TRUE
```
