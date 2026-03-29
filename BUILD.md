# GPBLUP Build Instructions

## Quick Start

```bash
cd /home/dhlee/GPBLUP
make clean
make all
# Executables: /home/dhlee/GPBLUP/bin/popQC, /home/dhlee/GPBLUP/bin/ReadFR
```

---

## Detailed Build Procedure

### Prerequisites

- **Compiler**: gfortran (GNU Fortran 90/95 compatible)
- **Build Tool**: make or CMake
- **System**: Linux (tested on Ubuntu/Debian)
- **Disk Space**: ~50 MB for build artifacts

### Verify Installation

```bash
# Check gfortran availability
gfortran --version

# Output should show: GNU Fortran (GCC) X.X.X

# Check working directory
ls /home/dhlee/GPBLUP/src/M_Kinds.f90
```

---

## Build Architecture

### Module Compilation Chain

The project uses a two-tier compilation strategy:

```
Tier 1: Shared Modules (src/)
├─ M_Kinds.f90          [Step 1] - Data type definitions
└─ M_Stamp.f90          [Step 2] - Version/timestamp utilities
         ↓
       build/common/ (shared object files)
         ↓
Tier 2: Unit Programs
├─ popQC.f90            [Step 3a] + M_Kinds.o + M_Stamp.o → bin/popQC
└─ ReadFR.f90           [Step 3b] + M_Kinds.o + M_Stamp.o → bin/ReadFR
```

### Compilation Sequence

**Step 1: Compile M_Kinds module**
```bash
cd /home/dhlee/GPBLUP/build/common
gfortran -c ../../src/M_Kinds.f90 \
    -ffree-line-length-none \
    -J . -I . \
    -g -O2
# Produces: M_Kinds.o, M_Kinds.mod
```

**Step 2: Compile M_Stamp module**
```bash
cd /home/dhlee/GPBLUP/build/common
gfortran -c ../../src/M_Stamp.f90 \
    -ffree-line-length-none \
    -J . -I . \
    -g -O2
# Produces: M_Stamp.o, M_Stamp.mod
```

**Step 3a: Compile popQC program**
```bash
cd /home/dhlee/GPBLUP/build/popQC
gfortran -c ../../popQC/popQC.f90 \
    -ffree-line-length-none \
    -J . -I ../common \
    -g -O2
# Produces: popQC.o

# Link executable
gfortran -o ../../bin/popQC popQC.o ../common/M_Kinds.o ../common/M_Stamp.o
# Produces: /home/dhlee/GPBLUP/bin/popQC (115 KB)
```

**Step 3b: Compile ReadFR program** (or skip if not needed)
```bash
cd /home/dhlee/GPBLUP/build/ReadFR
gfortran -c ../../ReadFR/ReadFR.f90 \
    -ffree-line-length-none \
    -J . -I ../common \
    -g -O2
# Produces: ReadFR.o

# Link executable
gfortran -o ../../bin/ReadFR ReadFR.o ../common/M_Kinds.o ../common/M_Stamp.o
# Produces: /home/dhlee/GPBLUP/bin/ReadFR (86 KB)
```

---

## Using Makefile

### Default Targets

#### Build All
```bash
make -C /home/dhlee/GPBLUP all
```
Builds both popQC and ReadFR with all dependencies.

#### Build popQC Only
```bash
make -C /home/dhlee/GPBLUP popQC
```

#### Build ReadFR Only
```bash
make -C /home/dhlee/GPBLUP ReadFR
```

#### Clean Build Artifacts
```bash
make -C /home/dhlee/GPBLUP clean
```
Removes all `.o` and `.mod` files from `build/` directories.

#### Clean and Rebuild Everything
```bash
make -C /home/dhlee/GPBLUP clean
make -C /home/dhlee/GPBLUP all
```

#### Full Clean (including executables)
```bash
make -C /home/dhlee/GPBLUP distclean
```

### Makefile Variables

Override default settings:

```bash
# Use custom compiler
make -C /home/dhlee/GPBLUP FC=ifort popQC     # Intel Fortran (if available)

# Enable verbose output
make -C /home/dhlee/GPBLUP VERBOSE=1 all

# Custom optimization level
make -C /home/dhlee/GPBLUP OPTFLAGS="-O3 -march=native" popQC
```

---

## Using CMake

### Initial Configuration

```bash
cd /home/dhlee/GPBLUP/build
cmake .. -DCMAKE_BUILD_TYPE=Release
```

### Build with CMake

```bash
cmake --build . --config Release
# Or
make -C /home/dhlee/GPBLUP/build
```

### Clean CMake Build

```bash
rm -rf /home/dhlee/GPBLUP/build/CMakeFiles
cmake --build . --target clean
```

---

## Compiler Flags Explained

| Flag | Purpose | Notes |
|------|---------|-------|
| `-c` | Compile only (no linking) | Produces `.o` file |
| `-o` | Output file name | Specifies executable name |
| `-J <dir>` | Module output directory | Where `.mod` files go |
| `-I <dir>` | Module include directory | Where to find `.mod` files |
| `-ffree-line-length-none` | Allow long lines | Required for popQC.f90 |
| `-g` | Include debug symbols | Enables gdb debugging |
| `-O0` | No optimization | Faster compilation |
| `-O2` | Moderate optimization | Default recommended |
| `-O3` | Aggressive optimization | Slower compilation, faster runtime |
| `-Wall` | Show all warnings | Recommended for development |
| `-Werror` | Treat warnings as errors | Strict mode |

### Recommended Compilation Flags

**Development Build** (fast compilation):
```bash
gfortran -c file.f90 -J dir -I dir -g -O0 -Wall
```

**Production Build** (optimized):
```bash
gfortran -c file.f90 -J dir -I dir -g -O3 -march=native
```

**Debug Build** (with debugging):
```bash
gfortran -c file.f90 -J dir -I dir -g -O0 -fbacktrace -fbounds-check
```

---

## Troubleshooting

### Issue 1: "Module not found" Error

**Symptom**: 
```
Error: Can't find module file 'm_kinds.mod' for module 'm_kinds'
```

**Causes**:
1. M_Kinds.f90 not compiled first
2. Module interface not in include directory
3. Typo in module name (case-sensitive)

**Solution**:
```bash
# Clean and rebuild in correct order
make -C /home/dhlee/GPBLUP clean
make -C /home/dhlee/GPBLUP all

# Or manually with explicit paths
cd /home/dhlee/GPBLUP/build/common
gfortran -c ../../src/M_Kinds.f90 -ffree-line-length-none -J . -I .
```

### Issue 2: "Line too long" Error

**Symptom**:
```
Error: Line is too long
```

**Cause**: popQC.f90 contains lines >132 characters (Fortran default limit)

**Solution**: Add `-ffree-line-length-none` flag:
```bash
gfortran -c popQC.f90 -ffree-line-length-none -J . -I ../common
```

### Issue 3: Linking Error

**Symptom**:
```
undefined reference to 'm_kinds_'
```

**Causes**:
1. Object files (`.o`) not linked
2. Module object file not included in link command
3. Wrong object file path

**Solution**:
```bash
# Include both object and module files in link command
gfortran -o popQC popQC.o ../common/M_Kinds.o ../common/M_Stamp.o

# Verify object files exist
ls -la build/popQC/*.o
ls -la build/common/*.o
```

### Issue 4: Old Object Files Causing Problems

**Symptom**: Compilation "succeeds" but executable behaves incorrectly

**Cause**: Stale `.o` files from previous compilation with different code

**Solution**: Clean rebuild
```bash
make -C /home/dhlee/GPBLUP clean
make -C /home/dhlee/GPBLUP all
```

### Issue 5: Missing Source Files

**Symptom**:
```
No such file or directory: M_Kinds.f90
```

**Causes**:
1. Working directory wrong
2. Path typo
3. File was moved

**Solution**:
```bash
# Verify file locations
ls -la /home/dhlee/GPBLUP/src/M_Kinds.f90
ls -la /home/dhlee/GPBLUP/popQC/popQC.f90

# Build from root directory
cd /home/dhlee/GPBLUP
make all
```

---

## Verification

### Verify Executable Creation

```bash
# Check if executables were created
ls -lh /home/dhlee/GPBLUP/bin/popQC
ls -lh /home/dhlee/GPBLUP/bin/ReadFR

# Expected output:
# -rwxr-xr-x 1 dhlee dhlee 115K [date] /home/dhlee/GPBLUP/bin/popQC
# -rwxr-xr-x 1 dhlee dhlee  86K [date] /home/dhlee/GPBLUP/bin/ReadFR
```

### Test popQC Execution

```bash
# Test with no arguments (should show usage)
/home/dhlee/GPBLUP/bin/popQC

# Expected output:
# ***** popQC: Genomic QC Pipeline v1.0 *****
# Usage: popQC parameter_file.txt
# ...

# Test with parameter file
/home/dhlee/GPBLUP/bin/popQC /home/dhlee/GPBLUP/popQC/parameter_example.txt
```

### Check Build Artifacts

```bash
# List all object files
ls -la /home/dhlee/GPBLUP/build/common/*.o
ls -la /home/dhlee/GPBLUP/build/popQC/*.o

# List all module interfaces
ls -la /home/dhlee/GPBLUP/build/common/*.mod
```

---

## Performance Tuning

### Compiler Optimization Levels

| Level | Complexity | Speed | Compilation Time |
|-------|-----------|-------|------------------|
| `-O0` | None | Slow | <1 sec |
| `-O1` | Basic | Medium | 1-2 sec |
| `-O2` | Moderate | Fast | 2-5 sec |
| `-O3` | Aggressive | Very Fast | 10-30 sec |
| `-Ofast` | Unsafe | Fastest | 30-60 sec |

### Recommended for Production

```bash
# Build with production optimizations
gfortran -c popQC.f90 \
    -ffree-line-length-none \
    -O3 -march=native \
    -J build/popQC -I build/common

# Result: 5-10% faster execution, minimal compilation overhead
```

---

## Advanced Topics

### Cross-Compilation

Build for different architecture:
```bash
# Build for ARM (e.g., Raspberry Pi)
gfortran -c file.f90 -march=armv7-a
```

### Static Linking

Ensure all libraries are embedded:
```bash
gfortran -static -o program object_files.o
```

### Shared Library Generation

Create library for external programs:
```bash
# Compile modules as shared object
gfortran -shared -fPIC -c M_Kinds.f90

# Create shared library
gfortran -shared -fPIC *.o -o libdkblupf90.so
```

### Profile-Guided Optimization

For maximum performance:
```bash
# Step 1: Build with profiling
gfortran -fprofile-generate -O3 popQC.f90

# Step 2: Run with sample data
./popQC sample_data.txt

# Step 3: Rebuild with profile feedback
gfortran -fprofile-use -O3 popQC.f90
```

---

## Continuous Integration

### Automated Build Script

```bash
#!/bin/bash
# auto_build.sh

set -e  # Exit on error

echo "=== GPBLUP Automated Build ==="
cd /home/dhlee/GPBLUP

# Clean
make clean

# Build
echo "Building modules..."
make all

# Test
echo "Running basic tests..."
./bin/popQC > /dev/null && echo "✓ popQC executable OK" || echo "✗ popQC failed"
./bin/ReadFR > /dev/null && echo "✓ ReadFR executable OK" || echo "✗ ReadFR failed"

# Verify sizes
echo "Binary sizes:"
ls -lh bin/popQC bin/ReadFR | awk '{print $9, $5}'

echo "=== Build Complete ==="
```

### Git Pre-Commit Hook

```bash
#!/bin/bash
# .git/hooks/pre-commit

# Prevent commits if build fails
cd /home/dhlee/GPBLUP
make clean
if ! make all; then
    echo "Build failed. Commit aborted."
    exit 1
fi
```

---

## Directory Cleanup

### After Successful Build

```bash
# Keep only executables and source
# Safe to remove:
rm -rf /home/dhlee/GPBLUP/build/*/*.o
rm -rf /home/dhlee/GPBLUP/build/*/*.mod

# These are rebuilt automatically when needed
```

### Total Space Usage

| Component | Size |
|-----------|------|
| Source code (*.f90) | 150 KB |
| Build artifacts | 10 MB |
| Executables | 200 KB |
| Libraries | 500 KB |
| Documentation | 2 MB |
| **Total** | **~13 MB** |

---

## Maintenance Schedule

| Task | Frequency | Command |
|------|-----------|---------|
| Clean rebuild | Weekly | `make -C /home/dhlee/GPBLUP clean all` |
| Verify executables | Monthly | `ls -lh bin/` |
| Update documentation | Per release | Edit `.md` files, regenerate PDF |
| Archive build artifacts | Quarterly | `make -C /home/dhlee/GPBLUP distclean` |

---

## Support & Troubleshooting

### Get Help

1. Check this document for common issues
2. Review compiler output (error messages are usually descriptive)
3. Verify file paths and permissions
4. Try clean rebuild: `make clean && make all`

### Report Issues

If problems persist:
1. Document exact error message
2. List compiler version: `gfortran --version`
3. Show build command: `make VERBOSE=1 popQC`
4. Include relevant file excerpts

---

**Version**: 1.0  
**Last Updated**: 2024  
**Status**: Production Ready  
**Maintained By**: GPBLUP Development Team

