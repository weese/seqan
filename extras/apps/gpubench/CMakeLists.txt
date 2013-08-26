# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/gpubench/CMakeLists.txt
#
# CMakeLists.txt file for gpubench.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_gpubench)
message (STATUS "Configuring extras/apps/gpubench")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# Search for CUDA.
find_package (CUDA QUIET)

# Disable CUDA acceleration if CUDA was not found.
if (NOT CUDA_FOUND)
  add_definitions (-DCUDA_DISABLED=1)
  message (STATUS "  CUDA not found: not building gpubench.")
  return ()
endif (NOT CUDA_FOUND)

option (CUDA_DISABLED "Set to ON to disable CUDA acceleration." OFF)

if (CUDA_DISABLED)
  add_definitions (-DCUDA_DISABLED=1)
  message (STATUS "  CUDA disabled: not building gpubench.")
  return ()
endif (CUDA_DISABLED)

# ----------------------------------------------------------------------------
# CUDA Setup
# ----------------------------------------------------------------------------

# Activate debug symbols for device code and remove optimizations.
#list (APPEND CUDA_NVCC_FLAGS_DEBUG -G -O0)

# Build from Tesla architecture upwards.
list (APPEND CUDA_NVCC_FLAGS -arch sm_20)

# Produce an error on illegal calls, e.g. __host__ __device__ --> __host__.
#list (APPEND CUDA_NVCC_FLAGS -Xcudafe --diag_error -Xcudafe h_hd_illegal_call)

# Disable all warnings.
list (APPEND CUDA_NVCC_FLAGS --disable-warnings)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
if (CUDA_FOUND AND NOT CUDA_DISABLED)
  cuda_add_executable(gpubench gpubench.cu gpubench.h myers.h testset.h)
endif (CUDA_FOUND AND NOT CUDA_DISABLED)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (gpubench ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (gpubench)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install gpubench in ${PREFIX}/bin directory
install (TARGETS gpubench DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/gpubench for SeqAn release builds.
install (FILES LICENSE README DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (gpubench)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable gpubench in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} gpubench CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:gpubench")
  set (CPACK_PACKAGE_NAME "gpubench")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "GPU Bench")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (gpubench "gpubench")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:gpubench")