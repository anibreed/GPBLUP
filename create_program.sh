#!/bin/bash

# ============================================================================
# New Program Template Generator for GPBLUP Project
# 
# Usage:
#   ./create_program.sh program_name [description]
#
# Example:
#   ./create_program.sh MyProgram "Description of MyProgram"
# ============================================================================

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_error() {
    echo -e "${RED}❌${NC} $1"
    exit 1
}

print_success() {
    echo -e "${GREEN}✅${NC} $1"
}

print_info() {
    echo -e "${BLUE}ℹ️${NC} $1"
}

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 program_name [description]"
    echo ""
    echo "Example:"
    echo "  $0 MyProgram 'My new analysis program'"
    exit 1
fi

PROGRAM_NAME=$1
DESCRIPTION=${2:-"New GPBLUP program"}
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Validate program name
if ! [[ $PROGRAM_NAME =~ ^[A-Za-z][A-Za-z0-9]*$ ]]; then
    print_error "Program name must start with a letter and contain only alphanumeric characters"
fi

# Create program directory
PROGRAM_DIR="${PROJECT_DIR}/${PROGRAM_NAME}"
if [ -d "$PROGRAM_DIR" ]; then
    print_error "Directory ${PROGRAM_NAME} already exists"
fi

echo ""
echo "Creating new program: ${PROGRAM_NAME}"
echo "Description: ${DESCRIPTION}"
echo ""

# ============================================================================
# Create directory structure
# ============================================================================
mkdir -p "${PROGRAM_DIR}/test"
print_success "Created directory: ${PROGRAM_DIR}"

# ============================================================================
# Create main program file
# ============================================================================
cat > "${PROGRAM_DIR}/${PROGRAM_NAME}.f90" << 'PROGRAM_TEMPLATE'
program PROGRAM_PLACEHOLDER
  ! =========================================================================
  ! Program: PROGRAM_PLACEHOLDER
  ! 
  ! Description:
  !   DESCRIPTION_PLACEHOLDER
  !
  ! Author(s):
  ! Date:
  ! Version: 1.0
  !
  ! Modules used:
  !   - M_Kinds: Kind parameter definitions
  !   - M_Variables: Global variable types and definitions
  !   - M_readpar: Parameter file parser
  !
  ! =========================================================================
  
  use M_Kinds
  use M_Variables
  use M_readpar
  implicit none
  
  ! Program control parameters
  character(len=256) :: param_file
  integer :: n_args
  
  ! =========================================================================
  ! Program header
  ! =========================================================================
  call print_header()
  
  ! =========================================================================
  ! Command line arguments
  ! =========================================================================
  n_args = command_argument_count()
  
  if (n_args < 1) then
    call print_usage()
    stop
  end if
  
  call get_command_argument(1, param_file)
  
  ! =========================================================================
  ! Main program flow
  ! =========================================================================
  print *, ""
  print *, "Parameter file: ", trim(param_file)
  print *, ""
  
  ! 1. Load and parse parameter file
  call setup_parameters(param_file)
  
  ! 2. Main logic
  call main_process()
  
  ! 3. Generate output
  call generate_output()
  
  ! =========================================================================
  ! Completion
  ! =========================================================================
  print *, ""
  print *, "✅ Program completed successfully!"
  print *, ""
  
contains

  ! =========================================================================
  ! Subroutines
  ! =========================================================================
  
  subroutine print_header()
    print *, ""
    print *, "==========================================="
    print *, "  PROGRAM_PLACEHOLDER v1.0"
    print *, "==========================================="
    print *, ""
    print *, "DESCRIPTION_PLACEHOLDER"
    print *, ""
  end subroutine print_header
  
  subroutine print_usage()
    print *, ""
    print *, "Usage: PROGRAM_PLACEHOLDER parameter_file"
    print *, ""
    print *, "Example:"
    print *, "  ./PROGRAM_PLACEHOLDER parameter.txt"
    print *, ""
  end subroutine print_usage
  
  subroutine setup_parameters(param_file)
    ! Load and parse parameter file using M_readpar module
    character(len=*), intent(in) :: param_file
    logical :: file_exists
    
    inquire(file=trim(param_file), exist=file_exists)
    if (.not. file_exists) then
      print *, "ERROR: Parameter file not found: ", trim(param_file)
      stop
    end if
    
    ! Parse parameter file using M_readpar module
    ! This populates global FileInfo variables used by all programs
    call read_parameters(trim(param_file))
    
    print *, "✓ Parameters loaded successfully"
    
  end subroutine setup_parameters
  
  subroutine main_process()
    ! Main processing logic
    print *, "Starting main process..."
    
    ! TODO: Implement your program logic here
    ! You can use global variables from M_Variables:
    !   - PEDFile, MAPFile, SNPFile, GENOFile, DATAFile
    !   - OutputPrefix, QCThresholds
    
    print *, "✓ Main process completed"
    
  end subroutine main_process
  
  subroutine generate_output()
    ! Generate output and reports
    print *, "Generating output..."
    
    ! TODO: Generate output files using OutputPrefix
    
    print *, "✓ Output generated"
    
  end subroutine generate_output

end program PROGRAM_PLACEHOLDER
PROGRAM_TEMPLATE

# Replace placeholders
sed -i "s/PROGRAM_PLACEHOLDER/${PROGRAM_NAME}/g" "${PROGRAM_DIR}/${PROGRAM_NAME}.f90"
sed -i "s|DESCRIPTION_PLACEHOLDER|${DESCRIPTION}|g" "${PROGRAM_DIR}/${PROGRAM_NAME}.f90"

print_success "Created: ${PROGRAM_NAME}.f90"

# ============================================================================
# Create README
# ============================================================================
cat > "${PROGRAM_DIR}/README.md" << 'README_TEMPLATE'
# PROGRAM_PLACEHOLDER

## Description
DESCRIPTION_PLACEHOLDER

## Build
The program is built as part of the integrated GPBLUP project using CMake.

From the GPBLUP root directory:
```bash
make rebuild
# or
cmake -B build -DCMAKE_BUILD_TYPE=Release
cd build
make -j$(nproc)
```

The executable will be created in `bin/PROGRAM_PLACEHOLDER`

## Usage
```bash
./bin/PROGRAM_PLACEHOLDER parameter_file
```

## Parameter File Format
See the template in `test/parameter_example.txt`

## Modules Used
- **M_Kinds**: Kind parameter definitions
- **M_Variables**: Global variable types (FileInfo, PEDInfo, QCThresholds)
- **M_readpar**: Parameter file parser (PEDFILE, MAPFILE, etc.)

## Integration with GPBLUP
This program is integrated into the unified GPBLUP project and uses:
- Common source modules from `source/`
- Shared parameter parsing via M_readpar
- Standard file formats (PED, MAP, GENO files)

When modules in `source/` are modified, this program is automatically recompiled.

## Development
1. Modify `PROGRAM_PLACEHOLDER.f90`
2. Update parameter templates if new parameters are needed
3. Run `make` from GPBLUP root to rebuild
4. Test with data in `test/` directory
README_TEMPLATE

sed -i "s/PROGRAM_PLACEHOLDER/${PROGRAM_NAME}/g" "${PROGRAM_DIR}/README.md"
sed -i "s|DESCRIPTION_PLACEHOLDER|${DESCRIPTION}|g" "${PROGRAM_DIR}/README.md"

print_success "Created: README.md"

# ============================================================================
# Create parameter template
# ============================================================================
cat > "${PROGRAM_DIR}/test/parameter_example.txt" << 'PARAM_TEMPLATE'
# Parameter file for PROGRAM_PLACEHOLDER
# 
# This parameter file defines input files and analysis parameters
# for the PROGRAM_PLACEHOLDER program

# ============================================================================
# Input Files (required)
# ============================================================================

PEDFILE: pedigree.txt
  HEADER: 1
  DELIM: SPACE
  NO_VARIABLES: 8
  1 ANIMAL_ID, 2 BREED, 3 SIRE, 4 DAM, 5 SEX, 6 BDATE, 7 LOC, 8 ARN

MAPFILE: map.txt
  HEADER: 0
  DELIM: SPACE
  NO_VARIABLES: 4
  1 SNP_ID, 2 CHR, 3 POS, 4 ALLELE

# ============================================================================
# Optional Input Files
# ============================================================================

# GENOFILE: genotype.txt
#   HEADER: 1
#   DELIM: SPACE
#   NO_VARIABLES: 8
#   1 ANIMAL_ID, 2 BREED, 3 SIRE, 4 DAM, 5 SEX, 6 BDATE, 7 LOC, 8 GENO

# ============================================================================
# Output
# ============================================================================

OUTPUTPREFIX: output/PROGRAM_PLACEHOLDER

# ============================================================================
# Program-specific parameters
# ============================================================================
# Add your program-specific parameters here
# Examples:
# ANALYSIS_TYPE: QC
# THRESHOLD_VALUE: 0.95
# NUM_THREADS: 4
PARAM_TEMPLATE

sed -i "s/PROGRAM_PLACEHOLDER/${PROGRAM_NAME}/g" "${PROGRAM_DIR}/test/parameter_example.txt"

print_success "Created: test/parameter_example.txt"

# ============================================================================
# Display next steps
# ============================================================================
echo ""
echo -e "${GREEN}=========================================="
echo "  Program template created successfully! "
echo "==========================================${NC}"
echo ""
print_info "Next steps:"
echo "  1. Edit: ${PROGRAM_DIR}/${PROGRAM_NAME}.f90"
echo "     - Implement main_process() subroutine"
echo "     - Add program-specific logic"
echo ""
echo "  2. Update: ${PROGRAM_DIR}/README.md"
echo "     - Add detailed description"
echo "     - Document parameters and usage"
echo ""
echo "  3. Integration: Add to CMakeLists.txt (automatic)"
echo "     - Edit: CMakeLists.txt"
echo "     - Add section for ${PROGRAM_NAME} (see template)"
echo ""
echo "  4. Build:"
echo "     cd /home/dhlee/GPBLUP"
echo "     make rebuild"
echo ""
echo "  5. Run:"
echo "     ./bin/${PROGRAM_NAME} ${PROGRAM_DIR}/test/parameter_example.txt"
echo ""
echo -e "${GREEN}=========================================="
echo "  Program directory: ${PROGRAM_DIR}"
echo "==========================================${NC}"
echo ""
