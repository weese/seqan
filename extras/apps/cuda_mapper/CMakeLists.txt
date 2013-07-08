# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/cuda_mapper/CMakeLists.txt
#
# CMakeLists.txt file for cuda_mapper.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_cuda_mapper)
message (STATUS "Configuring extras/apps/cuda_mapper")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP)
find_package (SeqAn REQUIRED)

# Search for CUDA.
find_package (CUDA QUIET)

# Disable CUDA acceleration if CUDA was not found.
if (NOT CUDA_FOUND)
  add_definitions (-DCUDA_DISABLED=1)
  message (STATUS "  CUDA not found: building cuda_mapper without CUDA acceleration.")
endif (NOT CUDA_FOUND)

option (CUDA_DISABLED "Set to ON to disable CUDA acceleration." OFF)

if (CUDA_DISABLED)
  add_definitions (-DCUDA_DISABLED=1)
  message (STATUS "  CUDA disabled: building cuda_mapper without CUDA acceleration.")
endif (CUDA_DISABLED)

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  OpenMP not found: building cuda_mapper without multi-threading.")
endif (NOT OPENMP_FOUND)

# ----------------------------------------------------------------------------
# CUDA Setup
# ----------------------------------------------------------------------------

# Activate debug symbols for device code and remove optimizations.
list (APPEND CUDA_NVCC_FLAGS_DEBUG -G -O0 -lineinfo)

# Activate verbose flag in ptxas to show register usage.
list (APPEND CUDA_NVCC_FLAGS -Xptxas --verbose)

# Enforce nvcc to use at most 32 registers per thread.
#list (APPEND CUDA_NVCC_FLAGS --maxrregcount 32)

# Build from Kepler architecture upwards.
list (APPEND CUDA_NVCC_FLAGS -arch sm_35)

# Produce an error on illegal calls, e.g. __host__ __device__ --> __host__.
list (APPEND CUDA_NVCC_FLAGS -Xcudafe --diag_error -Xcudafe h_hd_illegal_call)

# Disable all warnings.
list (APPEND CUDA_NVCC_FLAGS --disable-warnings)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

if (CMAKE_COMPILER_IS_GNUCXX)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mpopcnt")
endif (CMAKE_COMPILER_IS_GNUCXX)

set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")

if (CMAKE_COMPILER_IS_GNUCXX AND APPLE)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mno-avx")
endif (CMAKE_COMPILER_IS_GNUCXX AND APPLE)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
# TODO(esiragusa): add indexer to install below, or better use masai indexer.
add_executable(cuda_indexer indexer.cpp
                            types.h
                            ../masai/options.h
                            genome.h
                            genome_index.h)

if (CUDA_FOUND AND NOT CUDA_DISABLED)
  cuda_add_executable(cuda_mapper mapper.cpp
                                  mapper.h
                                  mapper.cu
                                  mapper.cuh
                                  types.h
                                  ../masai/options.h
                                  ../masai/store/reads.h
                                  genome.h
                                  genome_index.h)
else (CUDA_FOUND AND NOT CUDA_DISABLED)
  add_executable(cuda_mapper mapper.cpp
                             mapper.h
                             types.h
                             ../masai/options.h
                             ../masai/store/reads.h
                             genome.h
                             genome_index.h)
endif (CUDA_FOUND AND NOT CUDA_DISABLED)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (cuda_indexer ${SEQAN_LIBRARIES})
target_link_libraries (cuda_mapper ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (cuda_mapper)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install cuda_mapper in ${PREFIX}/bin directory
install (TARGETS cuda_mapper
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/cuda_mapper for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (cuda_mapper)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable cuda_mapper in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} cuda_mapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")
  set (CPACK_PACKAGE_NAME "cuda_mapper")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "CUDA Mapper")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (cuda_mapper "cuda_mapper")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")
