#!/bin/bash
# GPBLUP synthetic-demo runner. Run from this directory:  ./run_demo.sh
set -e
cd "$(dirname "$0")"
GPBLUP="../bin/gpblup"
if [ ! -x "$GPBLUP" ]; then
  echo "gpblup binary not found at $GPBLUP — chmod +x ../bin/* first"; exit 1
fi
# (re)generate the synthetic data if numpy is available; otherwise use the shipped files
if python3 -c "import numpy" 2>/dev/null; then
  echo ">> Generating synthetic data (deterministic, seed 2026)..."
  python3 make_demo.py
fi
echo ">> Running GPBLUP single-step evaluation (700 animals, 480 genotyped, 800 SNP)..."
OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" "$GPBLUP" model_demo.par | tail -n 25
echo ">> Validating GEBV against the hidden true breeding values..."
python3 eval_demo.py
