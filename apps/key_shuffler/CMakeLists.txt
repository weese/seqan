# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/key_shuffler/CMakeLists.txt
#
# CMakeLists.txt file for key_shuffler.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_key_shuffler CXX)
message (STATUS "Configuring apps/key_shuffler")

set (SEQAN_APP_VERSION "0.1.0")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
endif ()
find_package(CryptoPP REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (${CryptoPP_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (key_shuffler key_shuffler.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (key_shuffler ${SEQAN_LIBRARIES})
target_link_libraries (key_shuffler ${CryptoPP_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Generate text and man page information with key_shuffler_doc target.
add_custom_command (OUTPUT README.key_shuffler.txt
                    COMMAND key_shuffler --export-help=txt > README.key_shuffler.txt
                    DEPENDS key_shuffler)
add_custom_command (OUTPUT key_shuffler.1
                    COMMAND key_shuffler --export-help=man > key_shuffler.1
                    DEPENDS key_shuffler)
add_custom_target (key_shuffler_doc DEPENDS README.key_shuffler.txt key_shuffler.1)

# Build documentation on installation.
install (CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build ${CMAKE_BINARY_DIR} --target key_shuffler_doc)")

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (key_shuffler)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install key_shuffler in ${PREFIX}/bin directory
install (TARGETS key_shuffler
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/key_shuffler for SeqAn release builds.
install (FILES LICENSE
               README
               ${CMAKE_CURRENT_BINARY_DIR}/README.key_shuffler.txt
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/key_shuffler.1
         DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------
	
# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} key_shuffler CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:key_shuffler")
  set (CPACK_PACKAGE_NAME "key_shuffler")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "key_shuffler - Swiss Army Knife")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (key_shuffler "key_shuffler")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:key_shuffler")

