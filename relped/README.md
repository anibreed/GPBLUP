# relped

Pedigree-based relationship matrix preparation.

## Purpose

- Renumber pedigree records to a valid parent-before-offspring order.
- Compute A (additive relationship) and A inverse matrices.

## Usage

```bash
relped <parameter_file>
```

See USER_GUIDE.md in the repository root for the pipeline context.

## Source Layout
- relped main and relped-only modules live under relped/ and relped/src.
- Shared core modules remain under src/ (used by multiple programs).
