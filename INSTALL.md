# GPBLUP Installation Guide

## Requirements

- Linux (tested on Ubuntu/Debian)
- gfortran (Fortran 90+)
- make, cmake
- BLAS/LAPACK (e.g., OpenBLAS)
- OpenMP runtime

## One-batch Install (Recommended)

Build and install all programs in one step.

```bash
cd /home/dhlee/GPBLUP
./install.sh
# or
./install.sh /opt/gpblup
```

Installed binaries:
- ReadFR
- popQC
- relped
- relgeno

By default, binaries are installed to `~/.local/bin`. Make sure it is in your PATH.

## Build Only (No Install)

```bash
cd /home/dhlee/GPBLUP
make -C popQC rebuild
make -C ReadFR rebuild
make -C relped rebuild
make -C relgeno rebuild
```

Binaries are created under `bin/`.

## Optional Environment Variables

- `RELPED_TRIPLET_EPS`: numeric threshold to skip tiny values in triplet outputs.
- `RELPED_LOG`: path to a log file for relped. If not set, no log file is created.

## Troubleshooting

- Missing BLAS/LAPACK: install `libopenblas-dev` or equivalent.
- gfortran not found: install `gfortran` via your package manager.
