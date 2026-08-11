# GPBLUP User Guide

**GPBLUP** is a matrix-free single-step genomic BLUP engine for multi-breed, multi-herd
pig genetic evaluation with genomic genetic groups. It is the terminal stage of the
genomQC → impQC → GPBLUP pipeline: it consumes curated genotypes and pedigree and
returns genomic estimated breeding values (GEBV), total merit, reliabilities, and
candidate rankings.

This distribution ships the **compiled executables**; the source code is maintained
privately and is not included.

---

## 1. Executables

| Binary | Role |
|--------|------|
| `bin/gpblup` | **Integrated engine** — runs the full evaluation from a single parameter file (recommended entry point). |
| `bin/renped` | Standalone pedigree renumber/recode utility. |
| `bin/gpblup_renum` | Renumber/recode stage (when run as separate steps). |
| `bin/gpblup_readpar` | Parameter-file parser / validator. |
| `bin/gpblup_blup` | Standalone BLUP solver stage. |

Most users only need `bin/gpblup`.

## 2. Running

```bash
chmod +x bin/*
./bin/gpblup  model.par
```

Threads are controlled by `OMP_NUM_THREADS`:

```bash
OMP_NUM_THREADS=16 ./bin/gpblup model.par
```

A complete, annotated example parameter file is in [`examples/model_demo.par`](../examples/model_demo.par).

## 3. Parameter file

The parameter file is a plain-text, section-based format:

```ini
[DATA]
PHENOTYPE:    phen.txt        # phenotype file
PHENO_HEADER: 1               # 1 = first line is a header
PHENO_DELIM:  SPACE           # SPACE | TAB | COMMA
PEDIGREE:     ped.txt         # pedigree file
GENOTYPE:     geno.txt        # genotype file

[TRAITS]
N_TRAIT: 2
OBS_COL: 5 6                  # phenotype columns holding the trait observations
MISSING: -999
NAMES:   ADG  BF

[MODEL]
# EFFECT = name  col  TYPE(FIXED|COV|RANDOM)  traits...  [STRUCT=...]
EFFECT = herd_year  3  FIXED   1 2
EFFECT = sex        4  FIXED   1 2
EFFECT = age        7  COV     1 2
EFFECT = animal     1  RANDOM  1 2  STRUCT=ANIMAL_H   # ANIMAL_H = single-step H kernel
EFFECT = litter     8  RANDOM  1    STRUCT=IID

[VARIANCES]
GENETIC  = 12.0 3.0 ; 3.0 9.0     # t×t (co)variance matrix, rows separated by ';'
RESIDUAL = 30.0 5.0 ; 5.0 25.0

[GENOMIC]
KERNEL:      H                 # H = single-step; A = pedigree only; G = genomic only
ENGINE:      APY               # APY (default) | DENSE | ssGTBLUP
APY_CORE:    4000              # APY core size; set 0 or 'AUTO' to derive from eigenvalue dimensionality
POLYGENIC_W: 0.05              # residual-polygenic proportion (tau)

[INFERENCE]
METHOD:   BLUP                 # BLUP | REML_AI
MAX_ITER: 100
CONV_TOL: 1e-10

[OUTPUT]
# EBV, total merit, reliabilities, candidate rankings, marker effects (toggled by keywords)
```

Notes:
- `APY_CORE: AUTO` sizes the core from the eigenvalue dimensionality of the genomic
  relationship matrix (the number of largest eigenvalues explaining a target fraction
  of the variance), matching the core to the genomic dimensionality of the population.
- `KERNEL: H` performs single-step GBLUP; the base-alignment scalars are estimated
  automatically from the data.
- For a multi-origin (multi-herd) population, enable the genomic genetic-group covariate
  so that animals of different commercial origin are placed on a common, unbiased scale
  and a total-merit output is produced.

## 4. Input file formats

**Pedigree** (`ped.txt`) — one animal per line: `animal  sire  dam` (unknown parents = 0).
Run `bin/renped` first, or let `bin/gpblup` renumber internally.

**Genotype** (`geno.txt`) — one animal per line: an animal ID followed by allele dosages
coded `0/1/2` (missing handled upstream by impQC; the panel is complete by construction).

**Phenotype** (`phen.txt`) — whitespace/TAB/comma-delimited columns; the parameter file
references columns by number (`OBS_COL`, and the `col` field of each `EFFECT`). An optional
header line is skipped when `PHENO_HEADER: 1`.

## 5. Outputs

Depending on the `[OUTPUT]` keywords, GPBLUP writes:
- estimated breeding values (GEBV) and, for multi-origin data, **total merit** (within-origin
  solution + origin contribution);
- individual **reliabilities** (Monte-Carlo, made consistent with the single-step predictor);
- **candidate rankings** and indirect predictions of genotyped selection candidates;
- **marker effects** with windowed variance summaries (weighted single-step GBLUP).

## 6. Pipeline context

GPBLUP consumes the curated genotypes and pedigree produced by the companion
[genomQC](https://github.com/anibreed/genomQC) (quality control and identity recovery) and
[impQC](https://github.com/anibreed/impQC) (imputation and reference-panel construction)
stages. Any equivalent inputs in the formats above can be substituted.
