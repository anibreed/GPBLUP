# PHimpute Parameter Guide

This document describes the PHimpute parameter file format and supported keys.

## Usage

```bash
phimpute <parameter_file>
```

## File Format

- Sections start with a header keyword ending in ':'
- Keys are case-insensitive.
- Comments start with # and are ignored.
- Paths can be quoted with single or double quotes.

Sections:
- PEDFILE:
- MAPFILE:
- CANDIDATE_GENOFILE:
- REFERENCE_GENOFILE: (optional)
- IMPUTE:
- MENDEL:
- HAPLOTYPE:
- LD:
- QC:
- OUTPUT:

## Global Options (outside sections)

```
STOP_ON_MISMATCH: 0|1
SPECIES: <preset name>
AUTOSOME_MAX: <int>
SEXCHR_X: <int>
SEXCHR_Y: <int>
SEXCHR_XY: <int>
SEXCHR_Z: <int>
SEXCHR_W: <int>
USE_ZW: 0|1
SEXCHR_OVERRIDE: AUTOSOME_MAX=18 X=23 Y=24 XY=25 Z=0 W=0 USE_ZW=0
```

QC options can also be set globally:

```
QC_ENABLE: 0|1
QC_SAMPLE_MISS_MAX: 0.10
QC_SNP_MISS_MAX: 0.02
QC_MIN_MAF: 0.01
QC_REMOVE_AMBIGUOUS: 0|1
QC_LOG_PREFIX: qc
```

## PEDFILE Section

Required fields: ANIMAL_ID (or ID), SIRE, DAM. Optional: SEX, BIRTH, BREED, LOCATION.

```
PEDFILE: /path/to/pedigree.txt
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 5
1 ANIMAL_ID
2 SIRE
3 DAM
4 SEX
5 BIRTH
```

## MAPFILE Section

Required fields: SNP_ID, CHR, POS.
Optional: ARRAY_ALL, ARRAY_CHR, ALLELE1, ALLELE2.

```
MAPFILE: /path/to/map.txt
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 5
1 SNP_ID
2 CHR
3 POS
4 ARRAY_ALL
5 ARRAY_CHR
```

## CANDIDATE_GENOFILE / REFERENCE_GENOFILE

Required fields: ANIMAL_ID, GENO.

```
CANDIDATE_GENOFILE: /path/to/candidate.geno
INPUT:
HEADER: 1
DELIM: SPACE
NO_VARIABLES: 2
1 ANIMAL_ID
2 GENO
```

REFERENCE_GENOFILE is optional and uses the same layout.

## IMPUTE Section

```
IMPUTE:
  MENDEL: 1
  PEDIGREE: 1
  HAPLOTYPE: 1
  LD: 1
```

Keys:
- MENDEL / FAMILY / FAMILY_PHASING
- PEDIGREE
- HAPLOTYPE (or HAP)
- LD / POPULATION

## MENDEL Section

```
MENDEL:
  ENABLE: 1
  PEDIGREE: 1
```

## HAPLOTYPE Section

```
HAPLOTYPE:
  ENABLE: 1
  TUNE: 1
  WINDOW: 24 28 32
  SHIFT: 8 10 12
```

## LD Section

```
LD:
  ENABLE: 1
  WINDOW: 2
  WINDOW_MAX: 5
```

## QC Section

```
QC:
  ENABLE: 0
  SAMPLE_MISS_MAX: 0.10
  SNP_MISS_MAX: 0.02
  MIN_MAF: 0.01
  REMOVE_AMBIGUOUS: 0
  LOG_PREFIX: qc
```

## OUTPUT Section

```
OUTPUT:
  GENO: /path/to/phimpute_imputed.geno
```

## Minimal Example

```
PEDFILE: PED_Total.txt
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 3
1 ANIMAL_ID
2 SIRE
3 DAM

MAPFILE: MAP_K.txt
INPUT:
HEADER: 1
DELIM: TAB
NO_VARIABLES: 3
1 SNP_ID
2 CHR
3 POS

CANDIDATE_GENOFILE: GENO_cleaned.geno
INPUT:
HEADER: 0
DELIM: SPACE
NO_VARIABLES: 2
1 ANIMAL_ID
2 GENO

IMPUTE:
  MENDEL: 1
  PEDIGREE: 1
  HAPLOTYPE: 1
  LD: 1

OUTPUT:
  GENO: phimpute_imputed.geno
```
