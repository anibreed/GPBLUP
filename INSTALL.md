# Installing GPBLUP

GPBLUP is distributed as **precompiled Linux (x86-64) executables**. There is nothing to build.

## 1. Requirements

| Component | Requirement |
|-----------|-------------|
| Operating system | Linux, x86-64 |
| Runtime libraries | glibc; OpenMP `libgomp`; LAPACK/BLAS (`liblapack`, `libblas`) |
| Python (optional) | Python 3 with Matplotlib — only for the plotting helper scripts |

Check the runtime libraries resolve:

```bash
ldd bin/gpblup     # should resolve libgomp, liblapack/libblas, libm, libc, ...
```

Install missing runtimes if needed:

```bash
# Debian/Ubuntu
sudo apt-get install libgomp1 liblapack3 libblas3
# RHEL/CentOS/Rocky
sudo dnf install libgomp lapack blas
```

## 2. Install

```bash
chmod +x bin/*
# optional: put the main engine on your PATH
sudo cp bin/gpblup /usr/local/bin/gpblup
```

## 3. Run

```bash
./bin/gpblup  examples/model_demo.par
OMP_NUM_THREADS=16 ./bin/gpblup examples/model_demo.par
```

See [docs/USER_GUIDE.md](docs/USER_GUIDE.md) for the parameter-file format, input/output
specifications, and the genomic genetic-group workflow.

## Binaries

| Binary | Role |
|--------|------|
| `gpblup` | Integrated engine (recommended entry point) |
| `renped` | Pedigree renumber/recode utility |
| `gpblup_renum` | Renumber/recode stage |
| `gpblup_readpar` | Parameter-file parser/validator |
| `gpblup_blup` | Standalone BLUP solver stage |

All binaries are the open, unprotected build (no node-locking or license enforcement).
