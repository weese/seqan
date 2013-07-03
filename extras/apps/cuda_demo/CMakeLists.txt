# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/cuda_demo/CMakeLists.txt
#
# CMakeLists.txt file for cuda_demo.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_cuda_demo)
message (STATUS "Configuring extras/apps/cuda_demo")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# Search for CUDA.
find_package (CUDA QUIET)

# Break out if CUDA was not found.
if (NOT CUDA_FOUND)
  message (STATUS "  CUDA not found: Not building cuda_demo.")
  return ()
endif (NOT CUDA_FOUND)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

option (FM_INDEX_32_BITS "Set to ON to use a 32 bits FM-index." OFF)

if (FM_INDEX_32_BITS)
    add_definitions (-DFM_INDEX_32_BITS=1)
endif (FM_INDEX_32_BITS)

# ----------------------------------------------------------------------------
# CUDA Setup
# ----------------------------------------------------------------------------

# Activate debug symbols for device code and remove optimizations.
list (APPEND CUDA_NVCC_FLAGS_DEBUG -G -O0)

# Build from Tesla architecture upwards.
list (APPEND CUDA_NVCC_FLAGS -arch sm_20)

# Produce an error on illegal calls, e.g. __host__ __device__ --> __host__.
list (APPEND CUDA_NVCC_FLAGS -Xcudafe --diag_error -Xcudafe h_hd_illegal_call)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
if (CUDA_FOUND)
  cuda_add_executable(cuda_demo cuda_demo.cpp cuda_demo.cu cuda_demo.h)
endif (CUDA_FOUND)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (cuda_demo ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (cuda_demo)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install cuda_demo in ${PREFIX}/bin directory
install (TARGETS cuda_demo
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/cuda_demo for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (cuda_demo)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable cuda_demo in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} cuda_demo CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_demo")
  set (CPACK_PACKAGE_NAME "cuda_demo")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "CUDA Finder")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (cuda_demo "cuda_demo")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_demo")
