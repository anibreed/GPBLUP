# GPBLUP User Guide

This guide covers the four main programs: ReadFR, popQC, relped, and relgeno.

## Common Notes

- All executables are placed under `bin/` (or `PREFIX/bin` if installed).
- Input files are controlled by parameter files described below.

## ReadFR

Purpose: QC and preprocessing for Illumina FinalReport files.

Usage:
```bash
ReadFR <parameter_file>
```

Key inputs typically include:
- FinalReport file
- Pedigree file
- Map file

## popQC

Purpose: SNP/animal QC pipeline.

Usage:
```bash
popQC <parameter_file>
```

## relped

Purpose: Pedigree renumbering and relationship matrices.

Usage:
```bash
relped <parameter_file>
```

Outputs:
- Renumbered pedigree file
- A and Ainv matrices (binary and/or text based on settings)

Notes:
- Text triplet files include `row col value original_row_id original_col_id`.
- Use `RELPED_TRIPLET_EPS` to skip tiny values in triplet outputs.
- Set `RELPED_LOG` to enable a run log file.

## relgeno

Purpose: Genomic relationship matrix processing.

Usage:
```bash
relgeno <parameter_file>
```

## Parameter Files

Each program has its own parameter file format. See the example files in the
project directories for exact field names and layouts.

## Output Locations

- Default: output paths are controlled by parameter files.
- Binaries: `bin/` or `PREFIX/bin` if installed.
