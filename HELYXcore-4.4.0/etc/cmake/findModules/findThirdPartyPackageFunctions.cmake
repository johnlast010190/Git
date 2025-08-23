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
    Functions used to find thirdParty packages

[----------------------------------------------------------------------------]]



# ============================================================================ #
# ------------------------------- Functions ---------------------------------- #
# ============================================================================ #

# Function to make a new ThirdParty library
# Sets private variables
function(initialise_ThirdParty_object library expected_libs expected_incs)
    # private (i.e. CMake variables for internal use)
    set(${library}_EXPECTED_LIBS    "${expected_libs}"  PARENT_SCOPE )
    set(${library}_FOUND_LIBS       ""                  PARENT_SCOPE ) # only libs name
    set(${library}_FOUND_LIBS_DIR   ""                  PARENT_SCOPE ) # libs path and name
    set(${library}_MISSING_LIBS     ""                  PARENT_SCOPE )

    set(${library}_EXPECTED_INCS    "${expected_incs}"  PARENT_SCOPE )
    set(${library}_FOUND_INCS       ""                  PARENT_SCOPE )
    set(${library}_MISSING_INCS     ""                  PARENT_SCOPE )

    # Results
    set(${library}_FOUND            FALSE           PARENT_SCOPE )
    set(${library}_FOUND_MESSAGE    "${library}"    PARENT_SCOPE )
endfunction()


function(find_third_party_libraries library)
    foreach(expected_library IN LISTS ${library}_EXPECTED_LIBS)
        find_library(temporary_library_path
            NAMES ${expected_library}
            HINTS ${${library}_ARCH_PATH}
            # Add MPI_INCLUDE_PATH to catch system MPI and system dependencies
            # thereof (e.g. ptscotch)
            PATHS ${MPI_INCLUDE_PATH}
            PATH_SUFFIXES lib lib/${HELYX_MPI_NAME} ${${library}_PATH_SUFFIXES}
            # The CMAKE_LIBRARY_PATH just includes the HELYX library output dir,
            # which we don't want to search
            NO_CMAKE_PATH
        )

        if(EXISTS "${temporary_library_path}")
            add_library(${expected_library} SHARED IMPORTED GLOBAL)
            if("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
                set_target_properties(${expected_library} PROPERTIES
                    IMPORTED_IMPLIB "${temporary_library_path}"
                )
            else()
                set_target_properties(${expected_library} PROPERTIES
                    IMPORTED_LOCATION "${temporary_library_path}"
                )
            endif()
            # Now ${library}_FOUND_LIBS is a list of target names that is used later on
            list(APPEND ${library}_FOUND_LIBS ${expected_library})
            set(${library}_FOUND_LIBS ${${library}_FOUND_LIBS} PARENT_SCOPE)
            # For better report
            list(APPEND ${library}_FOUND_LIBS_DIR ${temporary_library_path})
            set(${library}_FOUND_LIBS_DIR ${${library}_FOUND_LIBS_DIR} PARENT_SCOPE)
        else()
            list(APPEND ${library}_MISSING_LIBS ${expected_library})
            set(${library}_MISSING_LIBS ${${library}_MISSING_LIBS} PARENT_SCOPE)
        endif()
        unset(temporary_library_path CACHE)
    endforeach()
endfunction()


# Almost the same as the above file, but probably not worth making generic...
function(find_third_party_includes library)
    string(TOLOWER ${library} libraryLC)
    foreach(expected_include IN LISTS ${library}_EXPECTED_INCS)
        find_file(temporary_include_path
            NAMES ${expected_include}
            HINTS ${${library}_ARCH_PATH}
            PATHS ${MPI_LIBRARY_PATHS}  # Useful when using system MPI and system dependencies thereof (e.g. ptscotch)
            PATH_SUFFIXES inc include include/${libraryLC} include/${HELYX_MPI_NAME} ${${library}_PATH_SUFFIXES}
        )
        if(EXISTS "${temporary_include_path}")
            list(APPEND ${library}_FOUND_INCS ${temporary_include_path})
            set(${library}_FOUND_INCS ${${library}_FOUND_INCS} PARENT_SCOPE)
        else()
            list(APPEND ${library}_MISSING_INCS ${expected_include})
            set(${library}_MISSING_INCS ${${library}_MISSING_INCS} PARENT_SCOPE)
        endif()
        unset(temporary_include_path CACHE)
    endforeach()
endfunction()


function(set_library_found_message library)

    if("${${library}_MISSING_INCS}" STREQUAL "" AND
        "${${library}_MISSING_LIBS}" STREQUAL "")
        # Success message
        set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}")
        #foreach(lib IN LISTS ${library}_FOUND_LIBS)
        foreach(lib IN LISTS ${library}_FOUND_LIBS_DIR)
            set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}
    ${lib}")
        endforeach()

        set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}")
        foreach(inc IN LISTS ${library}_FOUND_INCS)
            set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}
    ${inc}")
        endforeach()
        set(${library}_FOUND TRUE PARENT_SCOPE )

    else()
        # Failure messages in this order:
        # Disabled by user (${library}_REQUIRED set to OFF)
        # ARCH_PATH not set
        # ARCH_PATH not valid
        # Failed to find some stuff on arch path
        if(NOT "${${library}_REQUIRED}" MATCHES "^ON$|^AUTO$")
            set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}  (${library}_REQUIRED set to ${${library}_REQUIRED})")

        # ARCH_PATH not set
        elseif(NOT DEFINED ${library}_ARCH_PATH OR "${${library}_ARCH_PATH}" STREQUAL "")
            set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}  (${library}_ARCH_PATH not set)")

        # ARCH_PATH isn't valid'
        elseif(NOT EXISTS "${${library}_ARCH_PATH}")
            set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}  (${library}_ARCH_PATH set to the following value, which is not a valid path: \"${${library}_ARCH_PATH}\")")

        # Missing includes and libs
        elseif(NOT "${${library}_MISSING_INCS}" STREQUAL "" OR NOT "${${library}_MISSING_LIBS}" STREQUAL "")

            # Missing includes
            if(NOT "${${library}_MISSING_INCS}" STREQUAL "")
                set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}\n    Did not find the required header file(s) ")
                foreach(inc IN LISTS ${library}_MISSING_INCS)
                    set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE} \"${inc}\",")
                endforeach()
                string(TOLOWER ${library} libraryLC)
                set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE} using the following path and suffixes:\n    ${${library}_ARCH_PATH}/<inc;include;include/${libraryLC};include/${HELYX_MPI_NAME};${${library}_PATH_SUFFIXES}>")
            endif()

            # Missing libraries
            if(NOT "${${library}_MISSING_LIBS}" STREQUAL "")
                set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}\n    Did not find the required libraries ")
                foreach(lib IN LISTS ${library}_MISSING_LIBS)
                    set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE} \"${lib}\",")
                endforeach()
                set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE} using the following path and suffixes:\n    ${${library}_ARCH_PATH}/<lib;lib/${HELYX_MPI_NAME};${${library}_PATH_SUFFIXES}>")
            endif()

        endif()  # End of message setting
    endif()  # End of success/failure selection
    set(${library}_FOUND_MESSAGE "${${library}_FOUND_MESSAGE}" PARENT_SCOPE)
endfunction()


# This function search for the expected libraries and headers for a given ThirdParty library,
# and set output messages regarding the third-parties. It also associate the expected headers,
# with the expected libraries, depracating the usage of THIRDPARTY_${library_name}_INC.
# Use above functions and expose some variables to the parent scope.
#
# N.B.  This functions also trim the found includes paths based on the value of "expected_headers".
# Some ThirParty packages drops theirs headers in some folder structure under include/, e.g.:
#   IMATH headers are in "**/include/imath/*h", and its own project expects the includes dir to be "**/include/imath"
#   TBB headers are in "**/include/tbb/*h", but its own project expects the includes dir to be "**/include"
#   ALEMBIC headers are in somewhere like "**/include/Alembic/**/*h" but its own project expects the includes dir to be "*/include"
# To deal with those differences and to make the found include paths to resemble what the ThirdParty library expects,
# the include dir of a library is given by removing the "expected_headers" string from the found header full path.
#
# IMATH:
#   full path to header: **/imath-3.1.1/include/imath/ImathMath.h
#   expected_headers: "ImathMath.h"
#   include path will be: **/imath-3.1.1/include/imath/
# TBB:
#   full path to header: **/oneTBB-2021.7.0/include/tbb/tbb.h
#   expected_headers: "tbb/tbb.h"
#   include path will be: **/oneTBB-2021.7.0/include/
# ALEMBIC:
#   full path to header: **/alembic-1.8.3/include/Alembic/Abc/All.h
#   expected_headers: "Alembic/Abc/All.h"
#   include path will be: **/alembic-1.8.3/include/
#
# p.s.  If we don't respect the ThirdParty library base include dir, we get include headers errors
# from the ThirdParty library code when building HELYX applications
#
# BUG:  'helyx_find_thirdparty_package' (and subsequent calling functions 'find_third_party_*')
# can find libraries/includes under the MPI_PATH. This is problematic if a system MPI
# (usually installed under /usr/) is used, because CMake can find the library in the
# the system standard paths (if available) instead of the ${library}_ARCH_PATH.
# Later, system paths can be fully removed by 'strip_system_paths_from_thirdparty_includes'
# but ${library}_FOUND_MESSAGE will not be updated.
#
# Possible solution:  Make sure that ${library}_ARCH_PATH has precedence over the MPI_*_PATH
# on 'find_third_party_*'
function(helyx_find_thirdparty_package library_name expected_libraries expected_headers)

    if(${ARGC} GREATER 3)
        message(FATAL_ERROR
        "helyx_find_thirdparty_package expects exactly three arguments, but received ${ARGC}
(Usually caused by attempting to pass unquoted lists as single arguments)
")
    endif()

    initialise_ThirdParty_object(${library_name} "${expected_libraries}" "${expected_headers}")
    if("${${library_name}_REQUIRED}" MATCHES "^ON$|^AUTO$")
        find_third_party_libraries(${library_name})
        find_third_party_includes(${library_name})
    elseif("${${library_name}_REQUIRED}" STREQUAL "OFF")
        set(${library_name}_MISSING_INCS ${library_name}_EXPECTED_INCS)
        set(${library_name}_MISSING_LIBS ${library_name}_EXPECTED_LIBS)
    else()
        # Shouldn't be possible to get here, but better safe than sorry...
        message(FATAL_ERROR "${library_name}_REQUIRED set to \"${${library_name}_REQUIRED}\"
        Accepted values are ON, OFF, or AUTO")
    endif()
    set_library_found_message(${library_name})

    set(${library_name}_FOUND         "${${library_name}_FOUND}"         PARENT_SCOPE)
    set(${library_name}_FOUND_MESSAGE "${${library_name}_FOUND_MESSAGE}" PARENT_SCOPE)

    # If found, expose libraries and include dirs to parent scope
    if(${library_name}_FOUND)
        set(THIRDPARTY_${library_name} "${${library_name}_FOUND_LIBS}" PARENT_SCOPE)
        set(THIRDPARTY_${library_name}_DIR "${${library_name}_FOUND_LIBS_DIR}" PARENT_SCOPE)
        foreach(individual_file IN LISTS ${library_name}_FOUND_INCS)
            #get_filename_component(directory ${individual_file} DIRECTORY)
            # Remove the file path sufix from the path.
            # e.g.: <path>/include/openvdb/openvdb.h -> <path>/include/
            # Most of the ThirdParties uses only one header, so this loop should be ok
            foreach(header IN LISTS expected_headers)
                if("${individual_file}" MATCHES "${header}")
                    string(REPLACE "${header}" "" directory ${individual_file})
                endif()
            endforeach()
            # Glue them onto all associated IMPORTED targets that were synthesised.
            # The INCS variables probably should stop being used in favour of just hooking up to the IMPORTED targets
            # as is normal practice.
            foreach(tgt IN LISTS ${library}_FOUND_LIBS)
                # Now it's an interface include, and we have an imported target, one can just do
                # `target_include_directories(my_thing PRIVATE scotch)` and you'll get the include directories and lib.
                # No more mucking about with separately having to attach the libraries or includes onto all users.
                # Trimming the large amounts of now-mostly-unnecessary code that does that is left as an exercise for
                # the reader ;).
                # The THIRDPARTY_BLAH variables still exist, but are just a list of target names as opposed to filenames.
                # The THIRDPARTY_BLAH_INCLUDES vars also still exist, but are probably never what you want since linking
                # to the imported target will also get you the headers. Most calls to helyx_additional_includes can
                # probably now go away.
                target_include_directories(${tgt} INTERFACE "${directory}")
            endforeach()

            list(APPEND THIRDPARTY_${library_name}_INC "${directory}")
        endforeach()
        list(REMOVE_DUPLICATES THIRDPARTY_${library_name}_INC)
        set(THIRDPARTY_${library_name}_INC "${THIRDPARTY_${library_name}_INC}" PARENT_SCOPE)
    else()
        # Probably not necessary, but better safe than sorry...
        unset(THIRDPARTY_${library_name}_INC    PARENT_SCOPE)
        unset(THIRDPARTY_${library_name}        PARENT_SCOPE)
    endif()
endfunction()


# Including system include paths can be problematic
# find_package() may return multiple include paths of a library if different versions are installed,
# leading ot header name clashes and therefore build crashing. Moreover, multiple versions of the same lib
# can also be found if they exist in the machine, which may lead to the wrong library being used in the linking time.
#
# To avoid header name clashes or wrong library version linking, we should remove system stardard paths from
# THIRDPARTY_${library}_INC and THIRDPARTY_${library}, and to make sure that only the correct paths
# set in the *_ARCH_PATH are being used.


# Gets the compiler and the linker standard paths of includes and libraries
function(get_system_standard_paths)

    # ps.: String replace "\n" with ";" is failing randomly. Therefore we first
    # replace "\n" with "," and at the end we replace "," with ";"

    # ---- Get compiler include paths ---- #

    set(compiler_include_paths "")
    set(stdErr "") # making sure stdErr is empty

    # TO DO: Consider Only calling this execute_process once - store results in a cache
    # variable.  If C/CXX compiler changes, then the cache will need refreshing
    # anyway.  Careful though - would have to explicitly re-run if the compiler
    # changes, because it's possible to change compiler without deleting old cache.
    execute_process (
        COMMAND "${CMAKE_CXX_COMPILER}" "-E" "-x" "c++" "-v" "/dev/null"
        ERROR_VARIABLE stdErr
        OUTPUT_VARIABLE stdOut
        )

    if(stdErr)
        # Get the paths
        set(compiler_include_absolute_paths "")
        set(pattern "include \<\.\.\.\> search starts here\:(.*)End of search list\.")
        string(REGEX MATCH ${pattern} compiler_include_paths ${stdErr})
        if (compiler_include_paths)
            string(REPLACE "\n " "," compiler_include_paths ${compiler_include_paths})
            string(REPLACE "\n" "," compiler_include_paths ${compiler_include_paths})
            string(REPLACE "," ";" compiler_include_paths ${compiler_include_paths})
            list(REMOVE_AT compiler_include_paths 0 -1) # because list is: "include <...> search starts here: ; ... ; End of search list."
            # Convert each relative path to the absolute path
            foreach(lib_directory ${compiler_include_paths})
                get_filename_component(directory "${lib_directory}" ABSOLUTE)
                list(APPEND compiler_include_absolute_paths "${directory}")
            endforeach()
        endif()
    endif()

    # ---- Get linker and compiler libraries paths ---- #

    set(libraries_search_relative_paths "")
    set(libraries_search_absolute_paths "")
    set(stdOut "") # making sure stdOut is empty

    execute_process (
        COMMAND "${CMAKE_CXX_COMPILER}" "-print-search-dirs" "/dev/null"
        ERROR_VARIABLE stdErr
        OUTPUT_VARIABLE stdOut
        )

    if(stdOut)
        #set(pattern "(?<=LIBRARY_PATH)(.*)(?=COLLECT_GCC_OPTIONS)") # cmake does not support lazy (?) quantifiers
        #set(pattern "LIBRARY_PATH=(.*)/[^\n]*/") # get everything starting with LIBRARY_PATH= until the first \n
        set(pattern "(libraries: =)(.*)[^\n]") # get everything starting with LIBRARY_PATH= until the first \n
        string(REGEX MATCH ${pattern} compiler_library_paths ${stdOut})
        if(compiler_library_paths)
            string(REPLACE "libraries: =" "" compiler_library_paths ${compiler_library_paths}) # because it starts with 'libraries: ='
            string(REPLACE ":" "," compiler_library_paths ${compiler_library_paths})
            string(REPLACE "\n" "," compiler_library_paths ${compiler_library_paths})
            string(REPLACE "," ";" compiler_library_paths ${compiler_library_paths})
            list(APPEND libraries_search_relative_paths "${compiler_library_paths}")
        endif()
    endif()

    # making sure stdOut and linker_search_paths are empty befone next evaluation
    set(stdOut "") # making sure stdOut is empty
    set(linker_search_paths "")

    # ld --verbose | grep SEARCH_DIR | tr -s \' \;\' \\\\012
    execute_process (
        COMMAND ld --verbose
        ERROR_VARIABLE stdErr
        OUTPUT_VARIABLE stdOut
        )

    if(stdOut)
        set(pattern "SEARCH_DIR(.*)SECTIONS" )
        string(REGEX MATCH ${pattern} linker_search_paths ${stdOut})
        if(linker_search_paths)
            # it is a list of 'SEARCH_DIR("=/usr/lib");'
            string(REPLACE "SEARCH_DIR(\"=" "" linker_search_paths ${linker_search_paths}) # need to remove 'SEARCH_DIR("='
            string(REPLACE "\")" "" linker_search_paths ${linker_search_paths})# meed to remove '")'
            string(REPLACE "\n" "," linker_search_paths ${linker_search_paths})
            string(REPLACE " " "," linker_search_paths ${linker_search_paths})
            string(REPLACE "," ";" linker_search_paths ${linker_search_paths})
            list(REMOVE_AT linker_search_paths -1) # because it ends with "SECTIONS"
            list(APPEND libraries_search_relative_paths "${linker_search_paths}")
        endif()
    endif()

    # making sure stdOut and linker_search_paths are empty befone next evaluation
    set(stdOut "") # making sure stdOut is empty
    set(linker_search_paths "")

    # ldconfig -v 2>/dev/null | grep -v ^$'\t'
    execute_process (
        #COMMAND ldconfig -v
        COMMAND bash -c "echo \$(ldconfig -v 2>/dev/null | grep -v \^\$\'\\t\')"
        ERROR_VARIABLE stdErr
        OUTPUT_VARIABLE stdOut
        )

    set(linker_search_paths ${stdOut})
    # Need to avoid replacements in and empty variable. This is a bit ugly but works.
    # Maybe we could use regex?
    if(linker_search_paths)
        # linker_search_paths can contain only new lines, which will brake remaining replaces
        string(REPLACE ":\n" "\n" linker_search_paths ${stdOut}) # it may end with ':\n"
        string(REPLACE "\n" "" linker_search_paths ${linker_search_paths}) # remove any new line
        if(linker_search_paths)
            # linker_search_paths can contain only whitespaces, which will brake remaining replaces
            string(REPLACE " " "" linker_search_paths ${linker_search_paths})
            if(linker_search_paths)
                string(REPLACE ":" ";" linker_search_paths ${linker_search_paths}) # this must be the last one
                # then append to libraries_search_relative_paths
                list(APPEND libraries_search_relative_paths "${linker_search_paths}")
            endif()
        endif()
    endif()

    # Convert each relative path to the absolute path
    foreach(lib_directory ${libraries_search_relative_paths})
        get_filename_component(directory "${lib_directory}" ABSOLUTE)
        list(APPEND libraries_search_absolute_paths "${directory}")
    endforeach()

    # Remove any duplicated path
    list(REMOVE_DUPLICATES libraries_search_absolute_paths)
    list(LENGTH libraries_search_absolute_paths size)

    # ---- Set the cache variables --- #

    # now set the cache variable
    set(LIBRARY_SEARCH_PATHS "${libraries_search_absolute_paths}" CACHE INTERNAL "")
    mark_as_advanced(LIBRARY_SEARCH_PATHS)
    set(INCLUDE_SEARCH_PATHS "${compiler_include_absolute_paths}" CACHE INTERNAL "")
    mark_as_advanced(INCLUDE_SEARCH_PATHS)

endfunction()


# Takes libs from THIRDPARTY_${library}, gets their paths,
# prepends $ORIGIN, and adds them to the CMAKE_INSTALL_PATH.
# Also removes system standard paths
#
# PS.: Although we expected this function to be obsolete for Optional ThirdParty packages
# because we are using IMPORTED targets. We still need to append the library
# paths relative to ORIGIN on RPATH.
# THIRDPARTY_${library} contains only libraries name, which results in empty paths.
# If we indeed need to append to rpath, then we should use THIRDPARTY_${library}_DIR
function(append_library_paths_to_rpaths library)

    if("" STREQUAL "${THIRDPARTY_${library}}")
        if("TRUE" STREQUAL "${MPI_AWARE_COMPILER}" AND "MPI" STREQUAL "${library}")
            string(CONCAT s
                "Not appending anything to the RUNPATH for \"${library}\", as "
                "the compiler is mpi-aware.\n"
                "Please ensure your MPI setup is correct at run-time."
            )
            message(CLEAN "${s}")
        else()
            message(WARNING "Not appending anything to the RUNPATH for \"${library}\", as THIRDPARTY_${library} is empty\n")
        endif()
    endif()

    # Including system paths in RPATH can be problematic, especially on systems
    # which don't use RUNPATH.
    # Should probably also look in /etc/ld.so.conf
    list(APPEND system_library_paths
        "/lib"
        "/lib64"
        "/usr/lib"
        "/usr/lib64"
        "/usr/local/lib64"
        "/usr/local/lib"
        "/usr/lib/x86_64-linux-gnu"  # Ubuntu uses this path
        "/usr/lib64/x86_64-linux-gnu"  # Ubuntu uses this path
    )

    # Add compiler and linker paths of libraries
    list(APPEND system_library_paths LIBRARY_SEARCH_PATHS)
    list(REMOVE_DUPLICATES system_library_paths)

    # TODO: We may need to remove the *_ARCH_PATH if it is set to any of the
    # default system paths (see fftw installed from the package manager)
    # This is not an issue for FLEX and ZLIB since they are not added to the RPATH

    foreach(individual_lib ${THIRDPARTY_${library}_DIR})

        get_filename_component(directory ${individual_lib} DIRECTORY)

        # Sometimes EVTK just specifies library names without paths.  Not
        # entirely sure why, but it doesn't seem to matter - the RPATH is still
        # set correctly.'
        if("" STREQUAL "${directory}")
            continue()
        endif()

        get_filename_component(directory "${directory}" ABSOLUTE)

        if(EXISTS "${directory}" AND NOT "${directory}" IN_LIST system_library_paths)
            unset(worrying_library CACHE)
            # Warn if libc++ or libstdc++ exist on RPATH
            find_library (worrying_library
                NAMES c++ stdc++
                HINTS "${directory}"
                NO_DEFAULT_PATH
                NO_CMAKE_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH
            )
            if(worrying_library)
                get_filename_component(lib_name ${individual_lib} NAME)
                string(CONCAT s
                    "In order to satisfy the \"${lib_name}\" requirement "
                    "at run-time, the following path was added to the RPATH:\n"
                    "\t\"${directory}\"\n"
                    "However, the following library was found on this path:\n"
                    "\t\"${worrying_library}\"\n"
                    "This may conflict with libraries that HELYX requires at "
                    "run-time.  We recommend that all dependencies are kept on "
                    "separate paths.\n"
                )
                message(WARNING "${s}")
            endif()
            unset(worrying_library CACHE)

            file(RELATIVE_PATH directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${directory})
            list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/${directory}")
            # some libs like libopenfoam, libparhipdecomp and libptscotchdecomp
            # are located in lib/<HELYX_MPI_NAME>/, which is one level down
            # CMAKE_RUNTIME_OUTPUT_DIRECTORY
            list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../${directory}")
            set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} PARENT_SCOPE)
        endif()
    endforeach()
endfunction()

# TODO:  Investigate something similar for system library paths - this doesn't
# quite work, because linker not invoked.  Might be a more elegant way to handle
# libraries?
function(strip_system_paths_from_thirdparty_includes library)

    list(APPEND system_include_paths
        "/include"
        "/usr/include"
        "/usr/local/include"
        "/usr/include/x86_64-linux-gnu"  # Ubuntu uses this path
    )

    # Add compiler and linker paths of libraries
    list(APPEND system_include_paths INCLUDE_SEARCH_PATHS)
    list(REMOVE_DUPLICATES system_include_paths)

    # Check againts: STORED_THIRDPARTY_${library}_INC
    # library paths are (always?) absolute
    set(STORED_THIRDPARTY_${library}_INC "${THIRDPARTY_${library}_INC}")
    unset(THIRDPARTY_${library}_INC)

    # TODO: Should I remove or not the system paths if ${library}_ARCH_PATH is empty?
    #       This is principally for FLEX and ZLIB because we usually don`t set the arch path for them
    #
    # TODO: To keep an *_ARCH_PATH that is set to be a system standard path can still lead to
    # conflics with other libraries, then the solution is to move the library outside any system paths

    set(REMOVED_THIRDPARTY_${library}_INC "")
    foreach(lib_path ${STORED_THIRDPARTY_${library}_INC})
        # Need to remove paths from THIRDPARTY_${library}_INC if they match system or compiler include paths,
        # except if the path is the ${library}_ARCH_PATH
        if(EXISTS "${lib_path}" AND "${lib_path}" IN_LIST system_include_paths)
            set(lib_arch_path "${${library}_ARCH_PATH}")
            if("" STREQUAL "${lib_arch_path}")
                #continue()
            elseif(NOT "" STREQUAL "${lib_arch_path}" AND lib_path MATCHES "(${lib_arch_path}/include)")
                #continue()
            else()
                list(APPEND REMOVED_THIRDPARTY_${library}_INC ${lib_path})
                list(REMOVE_ITEM STORED_THIRDPARTY_${library}_INC ${lib_path})
            endif()
            unset(lib_arch_path)
        endif()
    endforeach()

    # If THIRDPARTY_${library}_INC is empty, this can lead to compilation issues due to missing *.h files,
    # but this will probably happen only if the *_ARCH_PATH is not right.
    # THIRDPARTY_${library}_INC would be empty if the *_ARCH_PATH is wrong and the found library live in system paths.
    # Moreover, this will influence principally FLEX and ZLIB since they usually live on system paths
    # If we reach this point, additional warnings will be already displayed to the user
    string(CONCAT s
        "THIRDPARTY_${library}_INC is empty. \n"
        "${library} includes were found, but all of them are on system paths "
        "that are not in accordance with the ${library}_ARCH_PATH. "
        "This can lead to compilation issues due to missing header files. \n"
        "    Removed THIRDPARTY_${library}_INC:\n"
        "        \"${REMOVED_THIRDPARTY_${library}_INC}\" \n"
        "System paths were removed to avoid conflics between different libraries or versions. "
        "If these are the correct include paths, consider changing the ${library}_ARCH_PATH. \n"
        "We recommend that dependencies are kept on separate paths from the system standard path.\n"
    )
    if(${library}_FOUND AND "" STREQUAL "${STORED_THIRDPARTY_${library}_INC}")
        message(WARNING "${s}")
    endif()

    # OPTIONAL_THIRDPARTY_LIBRARIES variable is only defined after mandatory libraries (FLEX, ZLIB, ...) are called,
    # so we can use that to identify variables with different scopes...
    if(OPTIONAL_THIRDPARTY_LIBRARIES)
        # Optional libraries and includes are already exposed to parent scope
        set(THIRDPARTY_${library}_INC "${STORED_THIRDPARTY_${library}_INC}")
    else()
        # Mandatory libraries and includes need to be exposed to parent scope
        set(THIRDPARTY_${library}_INC "${STORED_THIRDPARTY_${library}_INC}" PARENT_SCOPE)
    endif()

endfunction()


# Function to compare a given path to the <library>_ARCH_PATH
# To be used as replacement of "if(given_path MATCHES <library>_ARCH_PATH)"
# because MATCHES is regex based
#
# BUG?
# Quoting path_to_compare (i.e "${path_to_compare}") converts all the content
# to one single string. If path_to_compare is a list of paths, basically only the first
# path of the list will be compared.
function(compare_to_arch_path library path_to_compare)
    set(MATCH_${library}_ARCH_PATH "" PARENT_SCOPE)
    string(LENGTH "${${library}_ARCH_PATH}" arch_path_lenght)
    string(SUBSTRING "${path_to_compare}" 0 "${arch_path_lenght}" path_to_compare_substring)
    string(COMPARE EQUAL "${path_to_compare_substring}" "${${library}_ARCH_PATH}" path_match)
    set(MATCH_${library}_ARCH_PATH "${path_match}" PARENT_SCOPE)
endfunction()
