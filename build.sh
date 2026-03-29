#!/bin/bash

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
BUILD_DIR="$PROJECT_DIR/build"
BUILD_TYPE="${1:-Release}"

echo "=========================================="
echo "  GPBLUP Build Script"
echo "=========================================="
echo ""
echo "Project Dir: $PROJECT_DIR"
echo "Build Dir: $BUILD_DIR"
echo "Build Type: $BUILD_TYPE"
echo ""

# Check if source files exist
echo "Checking source files..."
if [ ! -f "$PROJECT_DIR/src/M_Kinds.f90" ]; then
  echo "ERROR: src/M_Kinds.f90 not found"
  exit 1
fi
echo "✓ Source files found"
echo ""

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Remove CMake cache
rm -f CMakeCache.txt
rm -rf CMakeFiles

# Configure with CMake
echo "Configuring with CMake..."
cmake "$PROJECT_DIR" \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_VERBOSE_MAKEFILE=OFF

if [ $? -ne 0 ]; then
  echo ""
  echo "ERROR: CMake configuration failed"
  exit 1
fi

echo ""
echo "Building..."
cmake --build . -- -j$(nproc)

if [ $? -ne 0 ]; then
  echo ""
  echo "ERROR: Build failed"
  exit 1
fi

echo ""
echo "=========================================="
echo "  Build Complete!"
echo "=========================================="
echo ""
echo "Executables created:"
ls -lh bin/
echo ""
echo "Run test:"
echo "  $BUILD_DIR/bin/test_readpar"
echo ""
echo "Run main program:"
echo "  $BUILD_DIR/bin/ReadFR /home/dhlee/DATA/FinalReport/param"
echo ""