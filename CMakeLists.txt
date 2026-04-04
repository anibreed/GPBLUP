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

# 소스 파일 - 공통 모듈 자동 수집
file(GLOB SOURCES /home/dhlee/GPBLUP/src/*.f90)
list(FILTER SOURCES EXCLUDE REGEX ".*/M_relped_io\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/libsearch\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/link_list\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/colleau_mod\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/colleau_text\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/diag_colleau\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/renped\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/read_amat\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/M_relped_service\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/Qsort4\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/libsearch\\.f90$")
list(FILTER SOURCES EXCLUDE REGEX ".*/M_Qsort\\.f90$")

# Main executable (ReadFR)
if(EXISTS "/home/dhlee/GPBLUP/ReadFR/ReadFR.f90")
  add_executable(ReadFR 
    /home/dhlee/GPBLUP/ReadFR/ReadFR.f90
    ${SOURCES}
  )
  set_target_properties(ReadFR PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/ReadFR
  )
  target_link_libraries(ReadFR PRIVATE lapack blas)
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
  target_link_libraries(popQC PRIVATE lapack blas)
endif()

# Main executable (relped)
if(EXISTS "/home/dhlee/GPBLUP/relped/relped.f90")
  add_executable(relped
    /home/dhlee/GPBLUP/src/M_param.f90
    /home/dhlee/GPBLUP/src/M_Kinds.f90
    /home/dhlee/GPBLUP/src/M_Variables.f90
    /home/dhlee/GPBLUP/src/M_Stamp.f90
    /home/dhlee/GPBLUP/src/M_Hashtable.f90
    /home/dhlee/GPBLUP/src/M_StrEdit.f90
    /home/dhlee/GPBLUP/src/M_ReadFile.f90
    /home/dhlee/GPBLUP/src/M_renped_config.f90
    /home/dhlee/GPBLUP/src/link_heap.f90
    /home/dhlee/GPBLUP/src/inbreed.f90
    /home/dhlee/GPBLUP/src/M_readpar_relped.f90
    /home/dhlee/GPBLUP/src/M_readpar_popqc.f90
    /home/dhlee/GPBLUP/src/M_readpar.f90
    /home/dhlee/GPBLUP/src/M_ped.f90
    /home/dhlee/GPBLUP/src/M_Sort.f90
    /home/dhlee/GPBLUP/src/M_relped_service.f90
    /home/dhlee/GPBLUP/src/link_list.f90
    /home/dhlee/GPBLUP/src/colleau_mod.f90
    /home/dhlee/GPBLUP/src/colleau_text.f90
    /home/dhlee/GPBLUP/src/renped.f90
    /home/dhlee/GPBLUP/relped/relped.f90
  )
  set_target_properties(relped PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/relped
  )
  target_include_directories(relped PRIVATE ${CMAKE_BINARY_DIR}/modules/relped)
  target_compile_options(relped PRIVATE -fopenmp -cpp)
  target_link_options(relped PRIVATE -fopenmp)
endif()

