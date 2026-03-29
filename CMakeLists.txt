cmake_minimum_required(VERSION 3.10)
project(GPBLUP Fortran)

# Compiler setup
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check -Wall -Wextra -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops -ffree-line-length-none")

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: /home/dhlee/GPBLUP")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

# 소스 파일 - 공통 모듈 자동 수집 (테스트 메인 제외)
file(GLOB SOURCES /home/dhlee/GPBLUP/src/*.f90)
list(FILTER SOURCES EXCLUDE REGEX ".*/test_readpar\\.f90$")

# Test executable
add_executable(test_readpar 
  /home/dhlee/GPBLUP/src/test_readpar.f90
  ${SOURCES}
)
set_target_properties(test_readpar PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/test_readpar
)

# Main executable (ReadFR)
if(EXISTS "/home/dhlee/GPBLUP/ReadFR/ReadFR.f90")
  add_executable(ReadFR 
    /home/dhlee/GPBLUP/ReadFR/ReadFR.f90
    ${SOURCES}
  )
  set_target_properties(ReadFR PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/ReadFR
  )
endif()

# Enable testing
enable_testing()
add_test(NAME test_readpar COMMAND test_readpar)
