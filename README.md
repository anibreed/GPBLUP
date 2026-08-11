# GPBLUP

**A matrix-free single-step genomic BLUP engine for multi-breed, multi-herd pig genetic evaluation with genomic genetic groups.**

GPBLUP is the terminal stage of the genomQC → impQC → GPBLUP pipeline. It takes curated genotypes and pedigree and returns genomic estimated breeding values (GEBV), total merit, reliabilities, and candidate rankings — driven by a single, self-configuring parameter file.

This repository distributes the **compiled executables and documentation** under the MIT License. The Fortran source code is not part of this distribution.

---

## Highlights

- **Matrix-free, iteration-on-data** multi-trait ssGBLUP — memory scales with the data, not with a stored coefficient matrix, so very large pedigrees stay tractable on a single machine.
- **APY** inverse genomic relationship with an **automatic core size** derived from the eigenvalue dimensionality of the genomic relationship matrix.
- **Reliabilities consistent with the single-step predictor** (genotyped animals correctly receive higher reliabilities).
- **Genomic genetic groups + total merit** — places animals of different commercial origin on a common, unbiased scale using an ancestry coordinate estimated directly from genotypes.
- **Optional weighted (non-linear) single-step GBLUP**; lossless pedigree reduction; indirect prediction of genotyped selection candidates.

## Requirements

- Linux (x86-64), OpenMP runtime (`libgomp`), LAPACK/BLAS runtime. See [INSTALL.md](INSTALL.md).

## Quick start

```bash
chmod +x bin/*
./bin/gpblup  examples/model_demo.par
```

`examples/model_demo.par` is an annotated example parameter file (two-trait single-step model, APY engine). See [docs/USER_GUIDE.md](docs/USER_GUIDE.md) for the parameter-file sections, input/output formats, and the genetic-group workflow.

> **Note:** a fully synthetic demonstration dataset (phenotype/pedigree/genotype files that run end-to-end with `model_demo.par`) will be added in a subsequent release. The example parameter file documents the required input formats.

## Data availability

The individual-level pig genotype and pedigree records used to validate GPBLUP are owned by a commercial breeding company and cannot be redistributed.

## Citation

See [CITATION.cff](CITATION.cff):

> Lee D. GPBLUP: a matrix-free single-step genomic BLUP engine for multi-breed, multi-herd pig genetic evaluation with genomic genetic groups. *Genetics Selection Evolution* (2026).
> Archived release: Zenodo, DOI 10.5281/zenodo.XXXXXXX

## License

MIT License — see [LICENSE](LICENSE). The compiled binaries and documentation are provided under the MIT License. The source code is maintained privately and is not part of this distribution.
