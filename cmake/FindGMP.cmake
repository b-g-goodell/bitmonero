# - Try to find libgmp
# Once done this will define
#
#  GMP_FOUND - system has libgmp
#  GMP_INCLUDE_DIR - the libgmp include directory
#  GMP_LIBRARIES - Link these to use libgmp
#  GMP_DEFINITIONS - Compiler switches required for using libgmp

# Copyright (c) 2006, Alexander Dymo, <adymo@kdevelop.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

find_path(GMP_INCLUDE_DIR gmp.h
  /usr/include
  /usr/local/include
)

find_library(GMP_LIBRARIES NAMES gmp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(gmp "Could not find libgmp" GMP_INCLUDE_DIR GMP_LIBRARIES)
# show the GMP_INCLUDE_DIR and GMP_LIBRARIES variables only in the advanced view
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES )

