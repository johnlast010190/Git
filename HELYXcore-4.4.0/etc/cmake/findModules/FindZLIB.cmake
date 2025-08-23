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
    (c) 2021 Engys Ltd.

Description
    Custom findModule for ZLIB

[----------------------------------------------------------------------------]]


# Push CMAKE_MODULE_PATH to enable calling CMake's find_package() from here.
set(STORED_CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "")

# First, we must make sure that ZLIB_ARCH_PATH is a real path for proper
# evaluations since it can assume relative paths.
if(NOT "" STREQUAL "${ZLIB_ARCH_PATH}")
    get_filename_component(ZLIB_ARCH_PATH ${ZLIB_ARCH_PATH} REALPATH)
endif()

# ps.: We need to force CMAKE_PREFIX_PATH to take only the ZLIB_ARCH_PATH for
# find_package, after that we can append everything to CMAKE_PREFIX_PATH.
# Need to do this because if system paths are defined for the ARCH_PATH of any
# other library (like system MPI in /usr/), that can take precedence and then
# ZLIB can be found outside ZLIB_ARCH_PATH.
# If ZLIB_ARCH_PATH is empty, then system path will be search anyway.

# Nothing special here - just respect ZLIB_ARCH_PATH
set(STORED_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
set(CMAKE_PREFIX_PATH "${ZLIB_ARCH_PATH}")

find_package(ZLIB QUIET REQUIRED)

set(CMAKE_PREFIX_PATH "${STORED_CMAKE_PREFIX_PATH}")
list(APPEND CMAKE_PREFIX_PATH "${ZLIB_ARCH_PATH}")

set(THIRDPARTY_ZLIB "${ZLIB_LIBRARIES}")
set(THIRDPARTY_ZLIB_DIR "${ZLIB_LIBRARIES}")
set(THIRDPARTY_ZLIB_INC "${ZLIB_INCLUDE_DIRS}")


# Set warning about ZLIB libs being found outside ZLIB_ARCH_PATH when
# ZLIB_ARCH_PATH is defined.
# Used when THIRDPARTY_ZLIB do not matches ZLIB_ARCH_PATH
string(CONCAT warn_message_libs
"ZLIB libraries found, but not on ZLIB_ARCH_PATH
    ZLIB_LIBRARIES:
        \"${THIRDPARTY_ZLIB_DIR}\"
    ZLIB_INCLUDES:
        \"${THIRDPARTY_ZLIB_INC}\"
    ZLIB_ARCH_PATH:
        \"${ZLIB_ARCH_PATH}\"
If these are the correct libraries, consider changing ZLIB_ARCH_PATH (or setting it to \"\")\n"
)

# Set warning about extra ZLIB being found
# when THIRDPARTY_ZLIB matches ZLIB_ARCH_PATH
string(CONCAT warn_message_extra_libs
"Multiple ZLIB libraries or includes were found in your system. \
This can leads to linker and runtime conflicts.
    ZLIB_LIBRARIES:
        \"${THIRDPARTY_ZLIB_DIR}\"
    ZLIB_INCLUDES:
        \"${THIRDPARTY_ZLIB_INC}\"
    ZLIB_ARCH_PATH:
        \"${ZLIB_ARCH_PATH}\"
If these are the correct libraries/includes, consider changing ZLIB_ARCH_PATH. \
We recommend that dependencies are kept on separate paths.\n"
)

list(LENGTH THIRDPARTY_ZLIB thirdParty_Zlib_size)
list(LENGTH THIRDPARTY_ZLIB_INC thirdParty_Zlib_inc_size)

if(NOT "" STREQUAL "${ZLIB_ARCH_PATH}")
    compare_to_arch_path("ZLIB" "${THIRDPARTY_ZLIB}")
    set(match_zlib_arch_path_lib "${MATCH_ZLIB_ARCH_PATH}")
    compare_to_arch_path("ZLIB" "${THIRDPARTY_ZLIB_INC}")
    set(match_zlib_arch_path_inc "${MATCH_ZLIB_ARCH_PATH}")
    if (NOT match_zlib_arch_path_lib OR
        NOT match_zlib_arch_path_inc)
        message(WARNING "${warn_message_libs}")
    endif()
else()
    if(thirdParty_zlib_size GREATER "1" OR thirdParty_zlib_inc_size GREATER "1")
        message(WARNING "${warn_message_extra_libs}")
    endif()
endif()

# Pop CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH "${STORED_CMAKE_MODULE_PATH}")
