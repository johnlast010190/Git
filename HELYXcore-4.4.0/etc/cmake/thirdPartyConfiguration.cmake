#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
[-----------------------------------------------------------------------------]
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2019-2021 Engys Ltd.

Description
    ThirdParty CMake configuration file for CMake

[----------------------------------------------------------------------------]]


message(TITLE "ThirdParty configuration")

include(${CMAKE_CURRENT_LIST_DIR}/findModules/findThirdPartyPackageFunctions.cmake)

get_system_standard_paths()

# ============================================================================ #
# ----------------------------- Required modules ----------------------------- #
# ============================================================================ #

# N.B.  Standard arguments to find_package() will propagate down the stack,
# so be careful!

if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    find_package(Regex)
    find_package(StackTrace)
endif ()

# Can't use MPI_REQUIRED, because the MPI FindPackageHandleStandardArgs sends a
# fatal error, which messes up our reporting.
find_package(MPI QUIET)
append_library_paths_to_rpaths(MPI)

# We must find Flex with find_package, because it supplies the flex_target macro
find_package(FLEX 2.5.35 QUIET REQUIRED)
strip_system_paths_from_thirdparty_includes(FLEX)

find_package(ZLIB QUIET REQUIRED)
strip_system_paths_from_thirdparty_includes(ZLIB)


# ============================================================================ #
# ----------------------------- Optional modules ----------------------------- #
# ============================================================================ #


include(${CMAKE_CURRENT_LIST_DIR}/parseThirdPartySettings.cmake)

foreach(library ${OPTIONAL_THIRDPARTY_LIBRARIES})
    initialise_search_variables(${library})
endforeach()


# ------------------------- Find optional libraries -------------------------- #

foreach(library ${OPTIONAL_THIRDPARTY_LIBRARIES})

    find_package(${library} QUIET)

    if(${${library}_FOUND})
        append_library_paths_to_rpaths(${library})
        strip_system_paths_from_thirdparty_includes(${library})
    endif()

endforeach()

# N.B.  Some libraries optionally link against the (pt)scotchDecomp library.
# We use a variable (${LIBRARY}DECOMP) to make that link optional.
#
# It looks like some libs linking against ${LIBRARY}DECOMP are being compiled before
# src/parallel/decompose, so we need ${LIBRARY}DECOMP to be defined here instead of
# in src/parallel/decompose/CMakeLists.txt
# Setting ${LIBRARY}DECOMP on src/parallel/decompose makes de decomposition methods
# to not be recognized by some applications as deomposePar

if(${SCOTCH_FOUND})
    set(SCOTCHDECOMP scotchDecomp)
endif()

if(${PTSCOTCH_FOUND})
    set(PTSCOTCHDECOMP ptscotchDecomp)
endif()

if(${KAHIP_FOUND})
    set(KAHIPDECOMP kahipDecomp)
endif()

if(${PARHIP_FOUND})
    set(PARHIPDECOMP parhipDecomp)
endif()

## This library is not currently supported
#
## helyx_find_thirdparty_package(BOOST "" "boost/config.hpp")


# ============================================================================ #
# -------------------------------- Messaging --------------------------------- #
# ============================================================================ #

message(CLEAN "")
message(SUBTITLE "Mandatory libraries")

# MPI is a special case: It's required, so we always assume it's found
message(CLEAN ${MPI_FOUND_MESSAGE})

# Same with Flex
message(CLEAN "FLEX
    ${THIRDPARTY_FLEX_DIR}
    ${THIRDPARTY_FLEX_INC}")

message(CLEAN "ZLIB
    ${THIRDPARTY_ZLIB_DIR}
    ${THIRDPARTY_ZLIB_INC}")

# And these libraries (when cross-compiling)
if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    message(CLEAN "REGEX
    ${THIRDPARTY_REGEX_DIR}
    ${THIRDPARTY_REGEX_INC}")
    message(CLEAN "STACK_TRACE
    ${THIRDPARTY_STACK_TRACE_DIR}
    ${THIRDPARTY_STACK_TRACE_INC}")
endif()

message(CLEAN "")


message(SUBTITLE "Optional libraries")
set(missing_libraries "")
foreach(library ${OPTIONAL_THIRDPARTY_LIBRARIES})
    if(${${library}_FOUND})
        message(CLEAN ${${library}_FOUND_MESSAGE})
    else()
        list(APPEND missing_libraries ${library})
    endif()
endforeach()

# HDF5 removed in version 3.4
thirdParty_package_deprecation_warning("HDF5" "3.4")

message(CLEAN "")

if(NOT "" STREQUAL "${missing_libraries}")
    message(SUBTITLE "Missing libraries")
    foreach(library ${missing_libraries})
        if("ON" STREQUAL ${library}_REQUIRED)
            set(MISSING_REQUIRED_LIB TRUE)
            message(SEND_ERROR "${${library}_FOUND_MESSAGE}")
        else()
            message(CLEAN "${${library}_FOUND_MESSAGE}")
        endif()
    endforeach()
    message(CLEAN "")
endif()
message(CLEAN
"[==============================================================================]

")




