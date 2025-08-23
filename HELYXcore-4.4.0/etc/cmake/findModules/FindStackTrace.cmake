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
    Custom findModule for StackTrace (used in the Windows build).

[----------------------------------------------------------------------------]]


# Must set arch paths to use the find_ThirdParty_package() function
# Only used for Windows, so ARCH_PATH can be hard-coded
set(STACK_TRACE_ARCH_PATH "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64GccDPInt32Opt/")

# ps.: Shouldn't use 'helyx_find_thirdparty_package' for STACK_TRACE as it is
# not in the OPTIONAL_THIRDPARTY_LIBRARIES list
#helyx_find_thirdparty_package(STACK_TRACE stack_trace stack_trace.h)



# ps.: Stack-Trace doesn't have a *Config.cmake file, so we can't call find_package()
# here as we do for ZLIB and FLEX. We need to manually find for the libs and
# incs, but we set the same message pattern.

# First, we must make sure that STACK_TRACE_ARCH_PATH is a real path for prope
# evaluations since it can assume relative paths.
if(NOT "" STREQUAL "${STACK_TRACE_ARCH_PATH}")
    get_filename_component(STACK_TRACE_ARCH_PATH ${STACK_TRACE_ARCH_PATH} REALPATH)
endif()

set(STACK_TRACE_FOUND FALSE)
set(THIRDPARTY_STACK_TRACE "")
set(THIRDPARTY_STACK_TRACE_INC "")

# find the library
find_library(temporary_library_path
    NAMES "stack_trace"
    HINTS "${STACK_TRACE_ARCH_PATH}"
    PATH_SUFFIXES lib
    # The CMAKE_LIBRARY_PATH just includes the HELYX library output dir,
    # which we don't want to search
    NO_CMAKE_PATH
)
# Import stack_trace and set important variables
if(EXISTS "${temporary_library_path}")
    add_library(stack_trace SHARED IMPORTED GLOBAL)
    set_target_properties(stack_trace PROPERTIES
        IMPORTED_IMPLIB "${temporary_library_path}"
    )
    set(THIRDPARTY_STACK_TRACE stack_trace) # contains only the lib name
    set(THIRDPARTY_STACK_TRACE_DIR "${temporary_library_path}") # contains the full path
    unset(temporary_library_path CACHE)
endif()

# Find the header
find_file(temporary_include_path
    NAMES "stack_trace.h"
    HINTS "${STACK_TRACE_ARCH_PATH}"
    PATH_SUFFIXES inc
)
# Remove the header file from the path
foreach(include_file IN LISTS temporary_include_path)
    get_filename_component(directory ${include_file} DIRECTORY)
    target_include_directories(stack_trace INTERFACE "${directory}")
    list(APPEND THIRDPARTY_STACK_TRACE_INC ${directory})
endforeach()
unset(temporary_include_path CACHE)

# Manually set STACK_TRACE_FOUND because we don't call find_package(), although it
# probably wont be used
if(${THIRDPARTY_STACK_TRACE} AND ${THIRDPARTY_STACK_TRACE_INC})
    set(STACK_TRACE_FOUND TRUE)
endif()



# Set warning about STACK_TRACE libs being found outside STACK_TRACE_ARCH_PATH
# when STACK_TRACE_ARCH_PATH is defined.
# Used when THIRDPARTY_STACK_TRACE do not matches STACK_TRACE_ARCH_PATH.
string(CONCAT warn_message_libs
"STACK_TRACE libraries found, but not on STACK_TRACE_ARCH_PATH
    STACK_TRACE_LIBRARIES:
        \"${THIRDPARTY_STACK_TRACE_DIR}\"
    STACK_TRACE_INCLUDES:
        \"${THIRDPARTY_STACK_TRACE_INC}\"
    STACK_TRACE_ARCH_PATH:
        \"${STACK_TRACE_ARCH_PATH}\"
If these are the correct libraries, consider changing STACK_TRACE_ARCH_PATH (or setting it to \"\")\n"
)

# Set warning about extra STACK_TRACE being found
# when THIRDPARTY_STACK_TRACE matches STACK_TRACE_ARCH_PATH
string(CONCAT warn_message_extra_libs
"Multiple STACK_TRACE libraries or includes were found in your system. \
This can leads to linker and runtime conflicts.
    STACK_TRACE_LIBRARIES:
        \"${THIRDPARTY_STACK_TRACE_DIR}\"
    STACK_TRACE_INCLUDES:
        \"${THIRDPARTY_STACK_TRACE_INC}\"
    STACK_TRACE_ARCH_PATH:
        \"${STACK_TRACE_ARCH_PATH}\"
If these are the correct libraries/includes, consider changing STACK_TRACE_ARCH_PATH. \
We recommend that dependencies are kept on separate paths.\n"
)

list(LENGTH THIRDPARTY_STACK_TRACE_DIR thirdParty_stacktrace_size)
list(LENGTH THIRDPARTY_STACK_TRACE_INC thirdParty_stacktrace_inc_size)

if(NOT "" STREQUAL "${STACK_TRACE_ARCH_PATH}")
    compare_to_arch_path("STACK_TRACE" "${THIRDPARTY_STACK_TRACE_DIR}")
    set(match_stack_trace_arch_path_lib "${MATCH_STACK_TRACE_ARCH_PATH}")
    compare_to_arch_path("STACK_TRACE" "${THIRDPARTY_STACK_TRACE_INC}")
    set(match_stack_trace_arch_path_inc "${MATCH_STACK_TRACE_ARCH_PATH}")
    if (NOT match_stack_trace_arch_path_lib OR
        NOT match_stack_trace_arch_path_inc)
        message(WARNING "${warn_message_libs}")
    endif()
else()
    if(thirdParty_stacktrace_size GREATER "1" OR
        thirdParty_stacktrace_inc_size GREATER "1")
        message(WARNING "${warn_message_extra_libs}")
    endif()
endif()
