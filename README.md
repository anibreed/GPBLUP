# GPBLUP

GPBLUP is a Fortran-based toolset for genomic data QC, pedigree processing, and
relationship matrix preparation. It provides multiple programs with shared
core modules and program-specific logic.

## Programs

- ReadFR: parse Illumina FinalReport files and emit GENO files
- popQC: SNP and animal QC pipeline
- relped: pedigree renumbering and A/Ainv matrix outputs
- relgeno: genomic relationship matrix processing
- phimpute: missing SNP imputation (pedigree-aware, configurable)

## Quick Start

Install all programs to a prefix:

```bash
./install.sh
# or
./install.sh /opt/gpblup
```

Build a single program from its directory:

```bash
make -C popQC rebuild
make -C ReadFR rebuild
make -C relped rebuild
make -C relgeno rebuild
```

Build everything with CMake:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Executables are placed in `bin/`.

## Documentation

- Build and install: [BUILD.md](BUILD.md), [INSTALL.md](INSTALL.md)
- Project layout: [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
- Program usage: [USER_GUIDE.md](USER_GUIDE.md)
- Algorithm background: [docs/ALGORITHM_BACKGROUND.md](docs/ALGORITHM_BACKGROUND.md)
- Model proof (sketches): [docs/MODEL_PROOF.md](docs/MODEL_PROOF.md)
- Packaging: [docs/PACKAGING.md](docs/PACKAGING.md)
- Full index: [DOCUMENTATION_INDEX.md](DOCUMENTATION_INDEX.md)

Program-specific docs:

- ReadFR: [ReadFR/README.md](ReadFR/README.md), [ReadFR/READFR_USER_MANUAL.md](ReadFR/READFR_USER_MANUAL.md)
- popQC: [popQC/README.md](popQC/README.md), [popQC/PIPELINE_GUIDE.md](popQC/PIPELINE_GUIDE.md)
- relped: [relped/README.md](relped/README.md)
- relgeno: [relgeno/RELGENO_PARAMS.md](relgeno/RELGENO_PARAMS.md)
- PHimpute: [PHimpute/README.md](PHimpute/README.md)

## Repository Layout (High Level)

```
src/        Shared modules
ReadFR/     ReadFR program and modules
popQC/      popQC program and modules
relped/     relped program and modules
relgeno/    relgeno program and modules
PHimpute/   phimpute program and modules
bin/        Executables
build/      Build artifacts
```

## Suggested Pipeline

1. ReadFR: FinalReport -> numeric GENO
2. popQC: population QC and filtering
3. phimpute: impute missing SNPs (pedigree-aware)
4. relgeno: G and G inverse
5. relped: A and A inverse