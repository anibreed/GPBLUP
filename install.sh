#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
PREFIX_ARG="${1:-}"

if [[ -n "$PREFIX_ARG" && "$PREFIX_ARG" != --* ]]; then
  PREFIX="$PREFIX_ARG"
else
  PREFIX="${PREFIX:-$HOME/.local}"
fi

BIN_DIR="$PROJECT_DIR/bin"
INSTALL_BIN="$PREFIX/bin"

echo "=========================================="
echo "  GPBLUP Install Script"
echo "=========================================="
echo "Project Dir: $PROJECT_DIR"
echo "Install Prefix: $PREFIX"
echo ""

mkdir -p "$BIN_DIR"

echo "Building popQC..."
make -C "$PROJECT_DIR/popQC" rebuild

echo "Building ReadFR..."
make -C "$PROJECT_DIR/ReadFR" rebuild

echo "Building relped..."
make -C "$PROJECT_DIR/relped" rebuild

echo "Building relgeno..."
make -C "$PROJECT_DIR/relgeno" rebuild

for exe in ReadFR popQC relped relgeno; do
  if [[ ! -x "$BIN_DIR/$exe" ]]; then
    echo "ERROR: Missing executable: $BIN_DIR/$exe"
    exit 1
  fi
done

mkdir -p "$INSTALL_BIN"
cp -f "$BIN_DIR/ReadFR" "$INSTALL_BIN/"
cp -f "$BIN_DIR/popQC" "$INSTALL_BIN/"
cp -f "$BIN_DIR/relped" "$INSTALL_BIN/"
cp -f "$BIN_DIR/relgeno" "$INSTALL_BIN/"

cat <<EOF

==========================================
  Install Complete!
==========================================
Installed binaries:
  $INSTALL_BIN/ReadFR
  $INSTALL_BIN/popQC
  $INSTALL_BIN/relped
  $INSTALL_BIN/relgeno
EOF
