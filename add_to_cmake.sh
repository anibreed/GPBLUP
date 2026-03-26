#!/bin/bash

# ============================================================================
# CMakeLists.txt Updater - Automatically add new program to build system
#
# Usage:
#   ./add_to_cmake.sh program_name
#
# Example:
#   ./add_to_cmake.sh MyProgram
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
    echo "Usage: $0 program_name"
    echo ""
    echo "Example:"
    echo "  $0 MyProgram"
    exit 1
fi

PROGRAM_NAME=$1
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CMAKE_FILE="${PROJECT_DIR}/CMakeLists.txt"

# Validate program name
if ! [[ $PROGRAM_NAME =~ ^[A-Za-z][A-Za-z0-9]*$ ]]; then
    print_error "Program name must start with a letter and contain only alphanumeric characters"
fi

# Check if program directory exists
if [ ! -d "${PROJECT_DIR}/${PROGRAM_NAME}" ]; then
    print_error "Directory ${PROGRAM_NAME} not found. Create it first with create_program.sh"
fi

# Check if source file exists
if [ ! -f "${PROJECT_DIR}/${PROGRAM_NAME}/${PROGRAM_NAME}.f90" ]; then
    print_error "Source file ${PROGRAM_NAME}.f90 not found"
fi

# Check if already in CMakeLists.txt
if grep -q "add_executable(${PROGRAM_NAME}" "${CMAKE_FILE}"; then
    print_error "Program ${PROGRAM_NAME} is already in CMakeLists.txt"
fi

echo ""
echo "Adding ${PROGRAM_NAME} to CMakeLists.txt..."
echo ""

# Create the CMake program block
read -r -d '' CMAKE_BLOCK << "EOF" || true

# ============================================================================
# PROGRAM_PLACEHOLDER
# ============================================================================
message(STATUS "Configuring PROGRAM_PLACEHOLDER program...")

set(PROGRAM_PLACEHOLDER_SOURCES
    PROGRAM_PLACEHOLDER/PROGRAM_PLACEHOLDER.f90
)

add_executable(PROGRAM_PLACEHOLDER ${PROGRAM_PLACEHOLDER_SOURCES})
target_link_libraries(PROGRAM_PLACEHOLDER PRIVATE gpblup_modules OpenMP::OpenMP_Fortran)
set_target_properties(PROGRAM_PLACEHOLDER PROPERTIES
    OUTPUT_NAME PROGRAM_PLACEHOLDER
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# PROGRAM_PLACEHOLDER depends on all modules
add_dependencies(PROGRAM_PLACEHOLDER gpblup_modules)

# Install to bin directory
install(TARGETS PROGRAM_PLACEHOLDER
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
)
EOF

# Replace placeholder with actual program name
CMAKE_BLOCK=${CMAKE_BLOCK//"PROGRAM_PLACEHOLDER"/"$PROGRAM_NAME"}

# Convert to uppercase for the variable name
PROGRAM_UPPER=$(echo "$PROGRAM_NAME" | tr '[:lower:]' '[:upper:]')
CMAKE_BLOCK=${CMAKE_BLOCK//"PROGRAM_PLACEHOLDER_SOURCES"/"${PROGRAM_UPPER}_SOURCES"}

# Find the last add_dependencies line and add after it
# But before the closing of the file
LAST_ADD_DEPS=$(grep -n "add_dependencies" "${CMAKE_FILE}" | tail -1 | cut -d: -f1)

if [ -z "$LAST_ADD_DEPS" ]; then
    print_error "Could not find add_dependencies in CMakeLists.txt"
fi

# Insert CMAKE_BLOCK after the last add_dependencies comment/message block
# Calculate the line to insert after (find next blank line or comment)
INSERT_LINE=$((LAST_ADD_DEPS + 30))  # Rough estimate

# Use sed to insert before the last "end module" or closing section
# Actually, let's append before the build_info section
BUILD_INFO_LINE=$(grep -n "custom target to display" "${CMAKE_FILE}" | head -1 | cut -d: -f1)

if [ -z "$BUILD_INFO_LINE" ]; then
    # Append before "add_custom_target"
    BUILD_INFO_LINE=$(grep -n "add_custom_target" "${CMAKE_FILE}" | head -1 | cut -d: -f1)
fi

if [ -z "$BUILD_INFO_LINE" ]; then
    print_error "Could not find insertion point in CMakeLists.txt"
fi

# Insert before the build_info section
INSERT_LINE=$((BUILD_INFO_LINE - 1))

# Create a temporary file
TEMP_FILE="${CMAKE_FILE}.tmp"

# Insert the new block
head -n "$INSERT_LINE" "${CMAKE_FILE}" > "$TEMP_FILE"
echo "" >> "$TEMP_FILE"
echo "$CMAKE_BLOCK" >> "$TEMP_FILE"
tail -n +$((INSERT_LINE + 1)) "${CMAKE_FILE}" >> "$TEMP_FILE"

# Replace original with temporary
mv "$TEMP_FILE" "${CMAKE_FILE}"

print_success "Added ${PROGRAM_NAME} to CMakeLists.txt"

# ============================================================================
# Verify and display next steps
# ============================================================================
echo ""
print_info "Verification:"

if grep -q "add_executable(${PROGRAM_NAME}" "${CMAKE_FILE}"; then
    echo "  ✓ add_executable found"
fi

if grep -q "add_dependencies(${PROGRAM_NAME}" "${CMAKE_FILE}"; then
    echo "  ✓ add_dependencies found"
fi

if grep -q "install(TARGETS ${PROGRAM_NAME}" "${CMAKE_FILE}"; then
    echo "  ✓ install target found"
fi

echo ""
print_info "Next steps:"
echo "  1. Build the new program:"
echo "     cd ${PROJECT_DIR}"
echo "     make rebuild"
echo ""
echo "  2. Verify binary was created:"
echo "     ls -lh bin/${PROGRAM_NAME}"
echo ""
echo "  3. Run the program:"
echo "     ./bin/${PROGRAM_NAME} ${PROGRAM_NAME}/test/parameter_example.txt"
echo ""

print_success "CMakeLists.txt update complete!"
