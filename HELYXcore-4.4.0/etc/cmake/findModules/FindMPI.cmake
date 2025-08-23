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
    (c) 2019-2020 Engys Ltd.

Description
    Custom findModule for MPI

[----------------------------------------------------------------------------]]


# Push CMAKE_MODULE_PATH to enable calling CMake's find_package() from here
set(STORED_CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "")

# First, we must make sure that MPI_ARCH_PATH is a real path for proper evaluations
# since it can assume relative paths.
if(NOT "" STREQUAL "${MPI_ARCH_PATH}")
    get_filename_component(MPI_ARCH_PATH ${MPI_ARCH_PATH} REALPATH)
endif()

# Then, we must explicitly ignore OPAL_PREFIX, which (if set incorrectly) can
# break MPI detection
unset(ENV{OPAL_PREFIX})
list(APPEND CMAKE_PREFIX_PATH "${MPI_ARCH_PATH}")

# To guess library version
set(MPI_DETERMINE_LIBRARY_VERSION TRUE CACHE INTERNAL "")

# For some reason, when compiling for Windows, CMake sometimes picks up system
# MPI libraries (e.g. libevent) and includes them in addition to the MS-MPI
# we're expecting from ThirdParty.  Therefore, don't use the CMake findModule
# when on Windows.
# It can also do this with headers
# This can be fixed (to a certain extent) with the following line, but I haven't
# done sufficient testing to re-instate the built-it findModule for the Windows
# build:
# set(MPI_GUESS_LIBRARY_NAME MSMPI CACHE INTERNAL "") # This also sets 'MPI_SKIP_COMPILER_WRAPPER' to 'true'
# Suspect this also helps:
#   ``MPI_ASSUME_NO_BUILTIN_MPI``
#       If true, the module assumes that the compiler itself does not provide an MPI
#       implementation and skips to step 2. Can be problematic at Fugaku
#   ``MPI_SKIP_COMPILER_WRAPPER``
#       If true, no compiler wrapper will be searched for. Can be problematic at Fugaku
#set(MPI_ASSUME_NO_BUILTIN_MPI  TRUE CACHE INTERNAL "")
#set(MPI_SKIP_COMPILER_WRAPPER  TRUE CACHE INTERNAL "")
#
# This could also help (Specify the base directory of the MPI installation):
# set(MPI_HOME "${MPI_ARCH_PATH}" CACHE INTERNAL "")
if(NOT "${HELYX_SYSTEM_NAME}" STREQUAL MSwindows)
    # Push and pop MPI_FIND_REQUIRED - we don't the built-in module to error
    set(STORED_MPI_FIND_REQUIRED "${MPI_FIND_REQUIRED}")
    set(MPI_FIND_REQUIRED "0")

    # Use QUIET to surpress any output from find_package since we already have
    # several warning/error messages related to MPI.
    # We can`t use find_package(MPI QUIET ${MPI_FIND_QUIETLY} HINTS "${MPI_ARCH_PATH}")
    # because it defaults to the config mode (it tried to find a MPIConfig.cmake file),
    # and then MPI can`t be find (we don`t use MPIConfig.cmake)
    find_package(MPI ${MPI_FIND_QUIETLY})

    # pop
    set(MPI_FIND_REQUIRED "${STORED_MPI_FIND_REQUIRED}")
endif()

if(${MPI_FOUND})
    # If MPI is found and MPI_CXX_COMPILER and CMAKE_CXX_COMPILER are the same
    # then, it is a mpi-aware compiler
    if("${MPI_CXX_COMPILER}" STREQUAL "${CMAKE_CXX_COMPILER}")
        set(MPI_AWARE_COMPILER TRUE)
    else()
        # The imported target made by the findmpi module.
        # The following properties are set among others:
        #   set_property(TARGET MPI::MPI_${LANG} PROPERTY INTERFACE_LINK_LIBRARIES "${MPI_${LANG}_LIBRARIES}")
        #   set_property(TARGET MPI::MPI_${LANG} PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_${LANG}_INCLUDE_DIRS}")
        set(THIRDPARTY_MPI MPI::MPI_CXX)

        # Set THIRDPARTY_MPI_DIR reporting only
        set(THIRDPARTY_MPI_DIR "${MPI_CXX_LIBRARIES}")

        # Set THIRDPARTY_MPI_INC for reporting and some includes
        # Different versions of CMake store the MPI headers in different places
        # I'm deliberately over-specifying THIRDPARTY_MPI_INC for robustness
        set(THIRDPARTY_MPI_INC ${MPI_C_INCLUDE_DIRS};${MPI_CXX_INCLUDE_DIRS};${MPI_C_HEADER_DIR};${MPI_CXX_HEADER_DIR};${MPI_INCLUDE_PATH};${MPI_C_INCLUDE_PATH};${MPI_CXX_INCLUDE_PATH})
        list(REMOVE_DUPLICATES THIRDPARTY_MPI_INC)

        # find_package can find multiple MPIs is they are available in the machine,
        # which can leat to wrong linking for the applications.
        # In this case, simply warn.
        list(LENGTH THIRDPARTY_MPI_DIR thirdparty_mpi_size)
        list(LENGTH THIRDPARTY_MPI_INC thirdparty_mpi_inc_size)
        if(thirdparty_mpi_size GREATER "1" OR thirdparty_mpi_inc_size GREATER "1")
            set(EXTRA_THIRDPARTY_MPI_DIR "")
            set(EXTRA_THIRDPARTY_MPI_INC "")
            foreach(lib_path ${MPI_CXX_LIBRARIES})
                compare_to_arch_path("MPI" "${lib_path}")
                if(NOT MATCH_MPI_ARCH_PATH)
                    list(APPEND EXTRA_THIRDPARTY_MPI_DIR ${lib_path})
                endif()
            endforeach()
            foreach(inc_path ${THIRDPARTY_MPI_INC})
                compare_to_arch_path("MPI" "${inc_path}")
                if(NOT MATCH_MPI_ARCH_PATH)
                    list(APPEND EXTRA_THIRDPARTY_MPI_INC ${inc_path})
                endif()
            endforeach()
        endif()
        if(EXTRA_THIRDPARTY_MPI_DIR)
            list(REMOVE_ITEM THIRDPARTY_MPI_DIR ${EXTRA_THIRDPARTY_MPI_DIR})
        endif()
        if(EXTRA_THIRDPARTY_MPI_INC)
            list(REMOVE_ITEM THIRDPARTY_MPI_INC ${EXTRA_THIRDPARTY_MPI_INC})
        endif()

        ## Rewrite MPI::MPI_CXX property to remove paths outside the MPI_ARCH_PATH
        #set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_LINK_LIBRARIES "${THIRDPARTY_MPI_DIR}")
        #set_property(TARGET MPI::MPI_CXX PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${THIRDPARTY_MPI_INC}")

        # What could possibly go wrong! (matches original semantics, though :D )
        target_include_directories(MPI::MPI_CXX INTERFACE ${THIRDPARTY_MPI_INC})
    endif()
else()
    set(MPI_REQUIRED ON)
    # Intel MPI keeps libmpi in a strange place, so this is necessary
    list(APPEND MPI_PATH_SUFFIXES lib/release)
    helyx_find_thirdparty_package(MPI mpi mpi.h)
endif()

# Set warning about MPI libs being found outside MPI_ARCH_PATH
# Used when THIRDPARTY_MPI_* do not matches MPI_ARCH_PATH
string(CONCAT warn_message_libs
"MPI libraries found, but not on MPI_ARCH_PATH.
    MPI_LIBRARIES:
        \"${THIRDPARTY_MPI_DIR}\"
    MPI_INCLUDES:
        \"${THIRDPARTY_MPI_INC}\"
    MPI_ARCH_PATH:
        \"${MPI_ARCH_PATH}\"
If these are the correct libraries, consider changing MPI_ARCH_PATH (or setting it to \"\")\n"
)

# Set warning about extra MPI libs being removed
# Used when THIRDPARTY_MPI matches MPI_ARCH_PATH but additional libs/includes were found
string(CONCAT warn_message_extra_libs
"MPI libraries/includes were found in addition to the MPI on MPI_ARCH_PATH.
    EXTRA MPI_LIBRARIES:
        \"${EXTRA_MPI_CXX_LIBRARIES}\"
    EXTRA MPI_INCLUDES:
        \"${EXTRA_THIRDPARTY_MPI_INC}\"
    MPI_ARCH_PATH:
        \"${MPI_ARCH_PATH}\"\n"
#"If these are the correct libraries, consider changing MPI_ARCH_PATH (or setting it to \"\")\n"
)

# for MPI-aware compilers, THIRDPARTY_MPI and THIRDPARTY_MPI_INC are empty
# so, find_package_handle_standard_args would return an error
if(NOT MPI_AWARE_COMPILER)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MPI
        REQUIRED_VARS THIRDPARTY_MPI THIRDPARTY_MPI_INC
        )
endif()

if(${MPI_FOUND})
    # Define MPI_FOUND_MESSAGE and set WARNING and ERROR messages
    if(MPI_AWARE_COMPILER)
        set(MPI_FOUND_MESSAGE "MPI\n    MPI-aware compiler: ${MPI_CXX_COMPILER}")
    else()
        if(MPI_CXX_LIBRARY_VERSION_STRING)
            string(REPLACE "," ";" MPI_CXX_LIBRARY_VERSION_LIST ${MPI_CXX_LIBRARY_VERSION_STRING})
            list(GET MPI_CXX_LIBRARY_VERSION_LIST 0 mpi_name_and_version)
            set(MPI_FOUND_MESSAGE "MPI: ${mpi_name_and_version}"
            "\n    ${THIRDPARTY_MPI_DIR}\n    ${THIRDPARTY_MPI_INC}")
        else()
            set(MPI_FOUND_MESSAGE "MPI\n    ${THIRDPARTY_MPI_DIR}\n    ${THIRDPARTY_MPI_INC}")
        endif()

        # At this stage, MPI must have been found

        # Set a list of directories that contain MPI libraries
        foreach(lib ${THIRDPARTY_MPI_DIR})
            get_filename_component(temp "${lib}" DIRECTORY)
            list(APPEND MPI_LIBRARY_PATHS "${temp}")
        endforeach()
        list(REMOVE_DUPLICATES MPI_LIBRARY_PATHS)

        # Check if mpirun or mpiexec are found
        if ("" STREQUAL "${MPIEXEC}")
            if (NOT "" STREQUAL "${MPIEXEC_EXECUTABLE}")
                set(MPIEXEC "${MPIEXEC_EXECUTABLE}")
            endif()

            find_program (MPIEXEC
                NAMES mpirun mpiexec
                HINTS ${MPI_ARCH_PATH} ${MPI_LIBRARY_PATHS}
                PATH_SUFFIXES bin;../bin
                DOC "Path on which mpiexec is found"
                )

            # We don't pack mpiexec.exe for windows. So, just warn for linux
            if (NOT EXISTS "${MPIEXEC}" AND NOT "${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
                message(WARNING
                    "Failed to find mpiexec, your PATH may be set incorrectly.  Please address this before running HELYX.")
            endif()
        endif()

        # Check if mpicc can build a simple program
        #check_mpicc("$MPI_ARCH_PATH}/bin/mpicc")

        # Extra warnings:
        # first check if MPI_ARCH_PATH is not in the list of found paths
        # then check if additional paths were found
        compare_to_arch_path("MPI" "${THIRDPARTY_MPI_DIR}")
        set(match_mpi_arch_path_lib "${MATCH_MPI_ARCH_PATH}")
        compare_to_arch_path("MPI" "${THIRDPARTY_MPI_INC}")
        set(match_mpi_arch_path_inc "${MATCH_MPI_ARCH_PATH}")
        if(NOT match_mpi_arch_path_lib OR NOT match_mpi_arch_path_inc)
            # only one path was found but it is not the arch path
            message(WARNING "${warn_message_libs}")
        elseif(EXTRA_THIRDPARTY_MPI_DIR OR EXTRA_THIRDPARTY_MPI_INC)
            # more paths were found in addition to the arch paths
            message(WARNING "${warn_message_extra_libs}")
        endif()

    endif()
else()
    set(MPI_FOUND_MESSAGE " MPI at\n\t${MPI_ARCH_PATH}")
    message(SEND_ERROR
        "Failed to find the following required library: \n${MPI_FOUND_MESSAGE}")
endif()

# Pop CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH "${STORED_CMAKE_MODULE_PATH}")
