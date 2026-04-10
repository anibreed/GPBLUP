#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIST_DIR="$ROOT_DIR/dist"
DATE_TAG="$(date +%Y%m%d)"
NAME="gpblup-${DATE_TAG}"
OUT_FILE="$DIST_DIR/${NAME}.tar.gz"

mkdir -p "$DIST_DIR"

EXCLUDES=(
  "--exclude=.git"
  "--exclude=dist"
  "--exclude=bin"
  "--exclude=build"
  "--exclude=build-*"
  "--exclude=*.o"
  "--exclude=*.mod"
  "--exclude=*.obj"
  "--exclude=*.a"
  "--exclude=*.so"
  "--exclude=*.dylib"
  "--exclude=*.exe"
  "--exclude=Geno_cleaned.txt"
  "--exclude=MAP_K.txt"
  "--exclude=MAP_K_sel.txt"
  "--exclude=PED_YY.txt"
)

tar -czf "$OUT_FILE" "${EXCLUDES[@]}" -C "$ROOT_DIR" .

echo "Wrote: $OUT_FILE"
