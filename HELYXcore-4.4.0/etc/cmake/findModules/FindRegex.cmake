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
    Custom findModule for REGEX (used in the Windows build).

[----------------------------------------------------------------------------]]


# Must set arch paths to use the find_ThirdParty_package() function
# Only used for Windows, so ARCH_PATH can be hard-coded
set(REGEX_ARCH_PATH "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64Gcc/mingw-libgnurx-2.5.1/")

# ps.: Shouldn't use 'helyx_find_thirdparty_package' for REGEX as it is not in
# the OPTIONAL_THIRDPARTY_LIBRARIES list
#helyx_find_thirdparty_package(REGEX regex regex.h)



# ps.: Regex doesn't have a *Config.cmake file, so we can't call find_package()
# here as we do for ZLIB and FLEX. We need to manually find for the libs and
# incs, but we set the same message pattern.

# First, we must make sure that REGEX_ARCH_PATH is a real path for proper
# evaluations since it can assume relative paths.
if(NOT "" STREQUAL "${REGEX_ARCH_PATH}")
    get_filename_component(REGEX_ARCH_PATH ${REGEX_ARCH_PATH} REALPATH)
endif()

set(REGEX_FOUND FALSE)
set(THIRDPARTY_REGEX "")
set(THIRDPARTY_REGEX_INC "")

# find the library
find_library(temporary_library_path
    NAMES "regex"
    HINTS "${REGEX_ARCH_PATH}"
    PATH_SUFFIXES lib
    # The CMAKE_LIBRARY_PATH just includes the HELYX library output dir,
    # which we don't want to search
    NO_CMAKE_PATH
)
# Import regex and set important variables
if(EXISTS "${temporary_library_path}")
    add_library(regex SHARED IMPORTED GLOBAL)
    set_target_properties(regex PROPERTIES
        IMPORTED_IMPLIB "${temporary_library_path}"
    )
    set(THIRDPARTY_REGEX regex) # contains only the lib name
    set(THIRDPARTY_REGEX_DIR "${temporary_library_path}") # contains the full path
    unset(temporary_library_path CACHE)
endif()

# Find the header
find_file(temporary_include_path
    NAMES "regex.h"
    HINTS "${REGEX_ARCH_PATH}"
    PATH_SUFFIXES include
)
# Remove the header file from the path and set the INTERFACE for the lib
foreach(include_file IN LISTS temporary_include_path)
    get_filename_component(directory ${include_file} DIRECTORY)
    target_include_directories(regex INTERFACE "${directory}")
    list(APPEND THIRDPARTY_REGEX_INC ${directory})
endforeach()
unset(temporary_include_path CACHE)

# Manually set REGEX_FOUND because we don't call find_package(), although it
# probably wont be used
if(${THIRDPARTY_REGEX} AND ${THIRDPARTY_REGEX_INC})
    set(REGEX_FOUND TRUE)
endif()



# Set warning about REGEX libs being found outside REGEX_ARCH_PATH when
# REGEX_ARCH_PATH is defined.
# Used when THIRDPARTY_REGEX do not matches REGEX_ARCH_PATH.
string(CONCAT warn_message_libs
"REGEX libraries found, but not on REGEX_ARCH_PATH
    REGEX_LIBRARIES:
        \"${THIRDPARTY_REGEX_DIR}\"
    REGEX_INCLUDES:
        \"${THIRDPARTY_REGEX_INC}\"
    REGEX_ARCH_PATH:
        \"${REGEX_ARCH_PATH}\"
If these are the correct libraries, consider changing REGEX_ARCH_PATH (or setting it to \"\")\n"
)

# Set warning about extra REGEX being found
# when THIRDPARTY_REGEX matches REGEX_ARCH_PATH
string(CONCAT warn_message_extra_libs
"Multiple REGEX libraries or includes were found in your system. \
This can leads to linker and runtime conflicts.
    REGEX_LIBRARIES:
        \"${THIRDPARTY_REGEX_DIR}\"
    REGEX_INCLUDES:
        \"${THIRDPARTY_REGEX_INC}\"
    REGEX_ARCH_PATH:
        \"${REGEX_ARCH_PATH}\"
If these are the correct libraries/includes, consider changing REGEX_ARCH_PATH. \
We recommend that dependencies are kept on separate paths.\n"
)

list(LENGTH THIRDPARTY_REGEX_DIR thirdParty_regex_size)
list(LENGTH THIRDPARTY_REGEX_INC thirdParty_regex_inc_size)

if(NOT "" STREQUAL "${REGEX_ARCH_PATH}")
    compare_to_arch_path("REGEX" "${THIRDPARTY_REGEX_DIR}")
    set(match_regex_arch_path_lib "${MATCH_REGEX_ARCH_PATH}")
    compare_to_arch_path("REGEX" "${THIRDPARTY_REGEX_INC}")
    set(match_regex_arch_path_inc "${MATCH_REGEX_ARCH_PATH}")
    if (NOT match_regex_arch_path_lib OR
        NOT match_regex_arch_path_inc)
        message(WARNING "${warn_message_libs}")
    endif()
else()
    if(thirdParty_regex_size GREATER "1" OR thirdParty_regex_inc_size GREATER "1")
        message(WARNING "${warn_message_extra_libs}")
    endif()
endif()
