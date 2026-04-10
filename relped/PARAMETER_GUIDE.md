# relped Parameter Guide

This document describes the relped parameter file format and supported keys.

## Usage

```bash
relped <parameter_file>
```

## File Format

- One setting per line.
- Keys are case-insensitive and typically end with a colon.
- Comments start with # and are ignored.
- Paths can be quoted with single or double quotes.

## Required Inputs

### PEDFILE

```
PEDFILE: /path/to/pedigree.txt
```

### PED Field Mapping

You can set the pedigree column indices in two ways.

Option A: Inline list

```
PEDINFO: 1 2 3 4 5
```

Order is:
1) ANIMAL_ID
2) SIRE_ID
3) DAM_ID
4) BIRTH (year/date)
5) SEX

Option B: INPUT block

```
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 5
1 ANIMAL_ID
2 SIRE_ID
3 DAM_ID
4 BIRTH
5 SEX
```

## Optional Reference Data

If you provide a reference file, relped can read additional ID/SEX info.

```
REFFILE: /path/to/reference.txt
REFINFO: 1 2
```

Order is:
1) ANIMAL_ID
2) SEX

Aliases:
- DATAFILE: same as REFFILE
- DATAINFO: same as REFINFO

## Outputs

### Renumbered pedigree file

```
FILENAME: /path/to/renumbered_ped.txt
```

### A and A inverse

```
RELATION: 1
RELATION_FILE: /path/to/A.mat
RELATION_INV: 1
RELATION_INV_FILE: /path/to/Ainv.mat
```

Output format:
- 1, BINARY, BIN (default)
- 2, TEXT, TXT, ASCII
- 3, BOTH (write both)

You can set format for both files:

```
OUTPUT_FORMAT: TEXT
```

Or per output block:

```
RELATIONSHIP: 1
FILE: /path/to/A.mat
FORMAT: TEXT

INV_RELATIONSHIP: 1
FILE: /path/to/Ainv.mat
FORMAT: TEXT
```

## Other Options

```
OMP_THREADS: 32
INBREED: 1
```

## Minimal Example

```
PEDFILE: PED_Total.txt
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 5
1 ANIMAL_ID
2 SIRE_ID
3 DAM_ID
4 BIRTH
5 SEX

RELATION: 1
RELATION_FILE: A.mat
RELATION_INV: 1
RELATION_INV_FILE: Ainv.mat
OUTPUT_FORMAT: TEXT

OMP_THREADS: 16
```
