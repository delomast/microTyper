cmake_minimum_required (VERSION 3.9)

# project name
project (microtyper)

find_package( ZLIB REQUIRED )
find_package(OpenMP)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# compiling included copy of bamtools (Version 2.5.1)
include(ExternalProject)

set(BAMTOOLS_ROOT ${CMAKE_CURRENT_BINARY_DIR}/bamtools)
set(BAMTOOLS_INCLUDE_DIRS "${BAMTOOLS_ROOT}/${CMAKE_INSTALL_INCLUDEDIR}/bamtools")
set(BAMTOOLS_LIBRARIES "${BAMTOOLS_ROOT}/${CMAKE_INSTALL_LIBDIR}")
set(bamtools_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/bamtools")
set(bamtools_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${bamtools_INSTALL_DIR})

include_directories(${BAMTOOLS_INCLUDE_DIRS})

ExternalProject_Add(bamtools
  PREFIX ${BAMTOOLS_ROOT}
  SOURCE_DIR ${BAMTOOLS_ROOT}
  BINARY_DIR ${BAMTOOLS_ROOT}/build
  INSTALL_DIR ${BAMTOOLS_ROOT}
  CMAKE_ARGS ${bamtools_CMAKE_ARGS}
)

add_library(libbamtools STATIC IMPORTED)
set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION ${BAMTOOLS_LIBRARIES}/libbamtools.a)
add_dependencies(libbamtools bamtools)

# likelihood calculations
add_executable(mtype src/microTyper.cpp src/utils.cpp src/utils.h src/calculate_llh.cpp src/calculate_llh_with_logSumExp.cpp src/calculate_llh.h src/typedefs.h src/count_matches.cpp)

add_executable(mtype2 src/microTyper_2.cpp src/utils.cpp src/utils.h src/countAlleleReads.cpp src/countAlleleReads.h src/typedefs.h)

include_directories(${BAMTOOLS_INCLUDE_DIRS})
target_link_libraries(mtype libbamtools ZLIB::ZLIB)
target_link_libraries(mtype2 libbamtools ZLIB::ZLIB)
# use openmp if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(mtype OpenMP::OpenMP_CXX)
    target_link_libraries(mtype2 OpenMP::OpenMP_CXX)
endif()

# genotype caller
add_executable(genoCaller src/llrGenotyper.cpp src/utils.cpp src/utils.h src/typedefs.h)

