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
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
seqan_add_cuda_executable(cuda_mapper mapper.cpp
                                      kernels.cu
                                      mapper.h
                                      kernels.h
                                      index.h
                                      ../masai/options.h
                                      ../masai/store/reads.h
                                      ../masai/store/genome.h
                                      ../masai/index/genome_index.h)


# Add dependencies found by find_package (SeqAn).
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
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")
  set (CPACK_PACKAGE_NAME "cuda_mapper")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "cuda_mapper")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")

  seqan_configure_cpack_app (cuda_mapper "cuda_mapper")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")

