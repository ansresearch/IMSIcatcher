# - Try to find fftw3 - the double-precision version of FFTW3
# Once done this will define
#  FFTW3_FOUND - System has fftw3
#  FFTW3_INCLUDE_DIRS - The fftw3 include directories
#  FFTW3_LIBRARIES - The libraries needed to use fftw3
#  FFTW3_DEFINITIONS - Compiler switches required for using fftw3

find_package(PkgConfig)
pkg_check_modules(PC_FFTW3 "fftw3 >= 3.0")
set(FFTW3_DEFINITIONS ${PC_FFTW3_CFLAGS_OTHER})

find_path(FFTW3_INCLUDE_DIR 
            NAMES fftw3.h
            HINTS ${PC_FFTW3_INCLUDEDIR} ${PC_FFTW3_INCLUDE_DIRS} $ENV{FFTW3_DIR}/include
            PATHS /usr/local/include 
                  /usr/include )

find_library(FFTW3_STATIC_LIBRARY
            NAMES fftw3.a libfftw3.a libfftw3-3.a
            HINTS ${PC_FFTW3_LIBDIR} ${PC_FFTW3_LIBRARY_DIRS} $ENV{FFTW3_DIR}/lib
            PATHS /usr/local/lib
                  /usr/lib)

find_library(FFTW3_LIBRARY 
            NAMES fftw3 libfftw3 libfftw3-3
            HINTS ${PC_FFTW3_LIBDIR} ${PC_FFTW3_LIBRARY_DIRS} $ENV{FFTW3_DIR}/lib
            PATHS /usr/local/lib
                  /usr/lib)

set(FFTW3_LIBRARIES ${FFTW3_LIBRARY} )
set(FFTW3_STATIC_LIBRARIES ${FFTW3_STATIC_LIBRARY} )
set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR} )

message(STATUS "FFTW3 LIBRARIES: " ${FFTW3_LIBRARIES})
message(STATUS "FFTW3 STATIC LIBRARIES: " ${FFTW3_STATIC_LIBRARIES})
message(STATUS "FFTW3 INCLUDE DIRS: " ${FFTW3_INCLUDE_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(fftw3  DEFAULT_MSG
                                  FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_STATIC_LIBRARY FFTW3_LIBRARY )
