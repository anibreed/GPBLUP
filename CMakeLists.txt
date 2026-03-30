cmake_minimum_required(VERSION 3.10)
project(GPBLUP Fortran)

# Compiler setup
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check -Wall -Wextra -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops -ffree-line-length-none")

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

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

# Main executable (popQC)
if(EXISTS "/home/dhlee/GPBLUP/popQC/popQC.f90")
  add_executable(popQC
    /home/dhlee/GPBLUP/popQC/popQC.f90
    ${SOURCES}
  )
  set_target_properties(popQC PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/popQC
  )
endif()

# Main executable (relped)
if(EXISTS "/home/dhlee/GPBLUP/relped/relped.f90")
  add_executable(relped
    /home/dhlee/GPBLUP/src/readpar_mod.f90
    /home/dhlee/GPBLUP/src/M_Kinds.f90
    /home/dhlee/GPBLUP/src/M_Variables.f90
    /home/dhlee/GPBLUP/src/M_Stamp.f90
    /home/dhlee/GPBLUP/src/M_Hashtable.f90
    /home/dhlee/GPBLUP/src/relped/M_param.f90
    /home/dhlee/GPBLUP/src/relped/M_ped.f90
    /home/dhlee/GPBLUP/src/relped/libsearch.f90
    /home/dhlee/GPBLUP/src/relped/M_readfile.f90
    /home/dhlee/GPBLUP/src/relped/link_list.f90
    /home/dhlee/GPBLUP/src/relped/colleau_mod.f90
    /home/dhlee/GPBLUP/src/relped/renped.f90
    /home/dhlee/GPBLUP/relped/relped.f90
  )
  set_target_properties(relped PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/relped
  )
  target_compile_options(relped PRIVATE -fopenmp -cpp)
  target_link_options(relped PRIVATE -fopenmp)
endif()

# Enable testing
enable_testing()
add_test(NAME test_readpar COMMAND test_readpar)
