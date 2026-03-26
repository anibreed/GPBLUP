cmake_minimum_required(VERSION 3.10)
project(GPBLUP Fortran)

# Compiler setup
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops")

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: /home/dhlee/GPBLUP")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

# 소스 파일 - 절대 경로 사용
set(SOURCES
  /home/dhlee/GPBLUP/source/M_kinds.f90
  /home/dhlee/GPBLUP/source/M_Variables.f90
  /home/dhlee/GPBLUP/source/M_StrEdit.f90
  /home/dhlee/GPBLUP/source/M_ReadFile.f90
  /home/dhlee/GPBLUP/source/M_readpar.f90
)

# Test executable
add_executable(test_readpar 
  /home/dhlee/GPBLUP/source/test_readpar.f90
  ${SOURCES}
)

# Main executable (ReadFR)
if(EXISTS "/home/dhlee/GPBLUP/source/ReadFR.f90")
  add_executable(ReadFR 
    /home/dhlee/GPBLUP/source/ReadFR.f90
    ${SOURCES}
  )
endif()

# Enable testing
enable_testing()
add_test(NAME test_readpar COMMAND test_readpar)
