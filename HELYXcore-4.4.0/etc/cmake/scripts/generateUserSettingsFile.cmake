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
    (c) 2019-2022 Engys Ltd.

Description
    This script generates a userSettings.cmake file from existing environment
    variables.

[----------------------------------------------------------------------------]]


# ============================================================================ #
# --------------------------------- Functions -------------------------------- #
# ============================================================================ #

function (normalise_path input_string output_variable)
    # Replace common directories (HELYX_PROJECT_DIR, WM_THIRD_PARTY_DIR, etc...)
    # with references to variables

    set(THIRDPARTY_PLATFORMS_PATH_LONG "$ENV{WM_THIRD_PARTY_DIR}/platforms/")
    string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG "linux")
    string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG "64")
    string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG "$ENV{WM_COMPILER}")
    set(THIRDPARTY_PLATFORMS_PATH_SHORT "${THIRDPARTY_PLATFORMS_PATH_LONG}")
    string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG "$ENV{WM_PRECISION_OPTION}")
    string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG "Int$ENV{WM_LABEL_SIZE}")

    set(normalised_input_dir "${input_string}")

    string(REPLACE "${THIRDPARTY_PLATFORMS_PATH_LONG}"
       "\${THIRDPARTY_PLATFORMS_PATH_LONG}" normalised_input_dir
       "${normalised_input_dir}"
       )
    string(REPLACE "${THIRDPARTY_PLATFORMS_PATH_SHORT}"
       "\${THIRDPARTY_PLATFORMS_PATH_SHORT}" normalised_input_dir
       "${normalised_input_dir}"
       )
    string(REPLACE "$ENV{WM_THIRD_PARTY_DIR}"
       "\${HELYX_THIRDPARTY_DIR}" normalised_input_dir
       "${normalised_input_dir}"
       )

    get_filename_component(normalised_project_dir "$ENV{WM_PROJECT_DIR}" ABSOLUTE)
    string(REPLACE "${normalised_project_dir}"
        "\${HELYX_PROJECT_DIR}" relative_path
        "${normalised_input_dir}"
    )

    get_filename_component(normalised_thirdParty_dir "$ENV{WM_THIRD_PARTY_DIR}" ABSOLUTE)
    string(REPLACE "${normalised_thirdParty_dir}"
        "\${HELYX_THIRDPARTY_DIR}" relative_path
        "${relative_path}"
    )

    # We could replace $ENV{HOME} as well, but not relying on environment
    # variables is preferable
    # get_filename_component(test "$ENV{HOME}" ABSOLUTE)
    # string(REPLACE "${test}"
    #     "\${HOME}" relative_path
    #     "${relative_path}"
    # )

    set(${output_variable} ${relative_path} PARENT_SCOPE)

endfunction()



# ============================================================================ #
# ---------------------------- Check WM settings ----------------------------- #
# ============================================================================ #

# Can't continue if the environment isn't set
set(BASHRC_ENVIRONMENT_SET TRUE)
# ThirdParty dir not strictly necessary
# CMAKE_C_COMPILER and CMAKE_CXX_COMPILER tested by cmake anyway
set(warning_string "")
foreach(expected_variable WM_PROJECT_DIR;WM_PROJECT_VERSION;WM_PRECISION_OPTION;WM_LABEL_SIZE;WM_COMPILE_OPTION)
    if(NOT DEFINED ENV{${expected_variable}})
        string(APPEND warning_string "${expected_variable}\n    ")
        set(BASHRC_ENVIRONMENT_SET FALSE)
    endif()
endforeach()

if("FALSE" STREQUAL "${BASHRC_ENVIRONMENT_SET}")
    message(FATAL_ERROR
        "
Environment not configured, userSettings.cmake will not be auto-generated
The following variables are required, and not set:
    ${warning_string}"
        )
endif()


# If the above variables are set, are they sensible?
# Warn here, but don't error.  Error later (in userSettingsConfiguration.cmake)
# if the CMake environment is wrong.
set(BASHRC_ENVIRONMENT_SENSIBLE TRUE)
set(warning_string "")

if(NOT IS_DIRECTORY "$ENV{WM_PROJECT_DIR}")
    string(APPEND warning_string
    "WM_PROJECT_DIR set to the following value, which is not a directory:
        $ENV{WM_PROJECT_DIR}
    "
    )
    set(BASHRC_ENVIRONMENT_SENSIBLE FALSE)
endif()

if(NOT "$ENV{WM_PRECISION_OPTION}" MATCHES "^SP|DP$")
    string(APPEND warning_string
    "WM_PRECISION_OPTION set to the following unrecognised value (should be SP or DP):
        $ENV{WM_PRECISION_OPTION}
    "
    )
    set(BASHRC_ENVIRONMENT_SENSIBLE FALSE)
endif()

if(NOT "$ENV{WM_LABEL_SIZE}" MATCHES "^32|64$")
    string(APPEND warning_string
    "WM_LABEL_SIZE set to the following unrecognised value (should be 32 or 64):
        $ENV{WM_LABEL_SIZE}
    "
    )
    set(BASHRC_ENVIRONMENT_SENSIBLE FALSE)
endif()

if(NOT "$ENV{WM_COMPILE_OPTION}" MATCHES "^Opt|Prof|Debug$")
    # Possibly not actually invalid - users *might* have implemented their
    # own compile options
    string(APPEND warning_string
    "WM_COMPILE_OPTION set to the following unrecognised value (should usually be Opt, Prof, or Debug):
        $ENV{WM_COMPILE_OPTION}
    "
    )
    set(BASHRC_ENVIRONMENT_SENSIBLE FALSE)
endif()

if("${BASHRC_ENVIRONMENT_SENSIBLE}" STREQUAL FALSE)
    message(WARNING "
Legacy environment set, but the following problems were found:
    ${warning_string}

    Attempting to generate userSettings.cmake anyway..."
    )
endif()



# ============================================================================ #
# ------------------------------ Generate file ------------------------------- #
# ============================================================================ #


# ------------------------------ Sanity checks ------------------------------- #

# script_location=`dirname $0`
# if [[ ! -d $ENV{WM_PROJECT_DIR} ]]; then
#     echo "Error:  \$ENV{WM_PROJECT_DIR} not a directory!
# Please source bashrc as if you were compiling with wmake before running this
# script"
#     exit 1
# elif [[ ! $ENV{WM_PROJECT_DIR} -ef $script_location/.. ]]; then
#     echo "Warning:  Expecting \$ENV{WM_PROJECT_DIR} equal to \
# `readlink -e $script_location/..`, but got the following value:
#     $ENV{WM_PROJECT_DIR}
# fi


# ----------------------------- File generation ------------------------------ #

set(output_file ${CMAKE_CURRENT_LIST_DIR}/../../userSettings.cmake)

file(WRITE ${output_file}
"#[[---------------------------------------------------------------------------]
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
    (c) 2020 Engys

Description
    This file encodes settings for compiling HELYXcore with CMake.

    This file can be automatically generated from the (now deprecated) wmake
    environment variables.  The wmake environment variables serve no other
    purpose, they are no longer used in the compilation of HELYXcore.

    The variables defined in this file follow the same conventions as the
    equivalent wmake environment variables or CMake build-in variables.  Where
    there is a conflict, the CMake style is preferred.

    Note that if \"FORCE\" is used for a variable defined in this file, the
    variable may no longer be set on the command line.  This ensures that this
    file remains the only place where variables may be set, and that the
    contents of this file will always be read in correctly when the CMake cache
    is refreshed.  Users who wish to set variables on the command line are
    welcome to remove the \"FORCE\" specification from this file, but should be
    aware that their cache may go out of date if not manually refreshed.

    If a required setting is not found in this file, a sensible default will be
    used.

[----------------------------------------------------------------------------]]



# ============================================================================ #
# ---------------------------- Directory Settings ---------------------------- #
# ============================================================================ #

set(HELYX_PROJECT_VERSION
    \"$ENV{WM_PROJECT_VERSION}\"
    CACHE STRING
    \"HELYX version\"
    FORCE
    )
set(HELYX_THIRDPARTY_VERSION
    \"$ENV{WM_THIRDPARTY_VERSION}\"
    CACHE STRING
    \"HELYX ThirdParty version\"
    FORCE
    )

# Get project path - this should work if this settings file is
# anywhere under etc/
get_filename_component(THIS_DIR \"\${CMAKE_CURRENT_LIST_DIR}\" REALPATH)
string(REGEX REPLACE \"^(.*)/etc.*$\" \"\\\\1/..\" PARENT_DIR \"\${THIS_DIR}\")
get_filename_component(PARENT_DIR \"\${PARENT_DIR}\" REALPATH)
set(HELYX_PROJECT_DIR
    \"\${PARENT_DIR}/HELYXcore-\${HELYX_PROJECT_VERSION}\"
    CACHE PATH
    \"Path to the top-level HELYX directory\"
    FORCE
    )
")


get_filename_component(TP_DIR "$ENV{WM_THIRD_PARTY_DIR}" REALPATH)
get_filename_component(TP_FROM_CORE_DIR
    "$ENV{WM_PROJECT_DIR}/../ThirdParty-$ENV{WM_THIRDPARTY_VERSION}" REALPATH
    )
if ("${TP_DIR}" STREQUAL "${TP_FROM_CORE_DIR}")
    set(TP_DIR_STRING "\"\${PARENT_DIR}/ThirdParty-\${HELYX_THIRDPARTY_VERSION}\"")
else()
    set(TP_DIR_STRING "\"$ENV{WM_THIRD_PARTY_DIR}\"")
endif()
file(APPEND ${output_file}
"set(HELYX_THIRDPARTY_DIR
    ${TP_DIR_STRING}
    CACHE PATH
    \"Path to the HELYX ThirdParty directory\"
    )
")

# Special treatment for ThirdParty compilers
set(gcc_string $ENV{WM_CC})
set(g++_string $ENV{WM_CXX})
foreach(compiler gcc g++)
    if ("$ENV{WM_COMPILER_TYPE}" STREQUAL "ThirdParty")
        find_program (third_party_compiler ${compiler}
            # Assume no-one's set the package root path variables, because this
            # switch new in CMake 3.12
            # NO_PACKAGE_ROOT_PATH

            # Don't search paths in cmake cache variables (e.g. CMAKE_PROGRAM_PATH)
            NO_CMAKE_PATH

            # Don't search paths in cmake env vars (e.g. $ENV{CMAKE_PROGRAM_PATH})
            NO_CMAKE_ENVIRONMENT_PATH

            # DO search on $ENV{PATH}
            # NO_SYSTEM_ENVIRONMENT_PATH

            # Don't search system paths
            NO_CMAKE_SYSTEM_PATH
            )
        if (NOT third_party_compiler)
            message(WARNING
                "\nWM_COMPILER_TYPE set to \"ThirdParty\", but CMake could not find the compiler \"$ENV{WM_CXX}\" (specified by WM_CXX) on the PATH environment variable.\nPlease check that the resulting user settings file contains the correct compiler information.
                ")
        else()
            get_filename_component(third_party_compiler "${third_party_compiler}" REALPATH)
            normalise_path(${third_party_compiler} ${compiler}_string)
            if ("${third_party_compiler}" STREQUAL "${${compiler}_string}")
                message(WARNING
                    "\nFailed to replace the following path with the HELYX_THIRDPARTY_DIR variable:\n\t${CXX_COMPILER_STRING}\nPlease check that the resulting user settings file contains the correct compiler information.
                ")
            endif()
        endif()
    endif()
endforeach()


if("$ENV{WM_COMPILER}" MATCHES "[A-Za-z]*\\d{1,4}")
    set(COMMENT "# ")
else()
    set(COMMENT "")
endif()

file(APPEND ${output_file}
"


# ============================================================================ #
# -----------------------------  Build Settings  ----------------------------- #
# ============================================================================ #

# If the executables specified by CMAKE_C_COMPILER and CMAKE_CXX_COMPILER are
# not available on the path, please specify the full paths to the executables.
set(CMAKE_C_COMPILER
    \"${gcc_string}\"
    CACHE PATH
    \"The C compiler to use (can be a path to an executable)\"
    FORCE
    )
set(CMAKE_CXX_COMPILER
     \"${g++_string}\"
    CACHE PATH
    \"The C++ compiler to use (can be a path to an executable)\"
    FORCE
    )

# HELYX_COMPILER_NAME is a user-friendly name for CMAKE_C_COMPILER.  It has no
# bearing on the compiler that will be used - it's only used to name
# subdirectories in the following locations:
#     \${HELYX_PROJECT_DIR}/platforms/,
#     \${HELYX_THIRDPARTY_DIR}/platforms/,
#     \${HELYX_PROJECT_DIR}/cbuild/.
# N.B.  For new builds, it is recommended that this variable is left unset.  If
# left unset, a suitable value will be automatically chosen.  However, it can be
# useful to set this variable to maintain backwards compatibility.
# HELYX_COMPILER_NAME is usually the same as CMAKE_C_COMPILER but with the first
# letter capitalised, sometimes followed by the first two digits of the compiler
# version number (e.g. \"Gcc\", \"Gcc71\", etc...)
# By default CMake will derive HELYX_COMPILER_NAME from CMAKE_C_COMPILER, but
# HELYX_COMPILER can be set manually as follows:
${COMMENT}set(HELYX_COMPILER_NAME
${COMMENT}    \"$ENV{WM_COMPILER}\"
${COMMENT}    CACHE PATH
${COMMENT}    \"User-friendly name for CMAKE_C_COMPILER (ONLY used for setting directory names)\"
${COMMENT}    FORCE
${COMMENT}    )

set(HELYX_PRECISION_OPTION
    \"$ENV{WM_PRECISION_OPTION}\"
    CACHE STRING
    \"Float precision, either double (\\\"DP\\\") or single precision (\\\"SP\\\")\"
    FORCE
    )
set(HELYX_LABEL_SIZE
    \"$ENV{WM_LABEL_SIZE}\"
    CACHE STRING
    \"Number of bits in an OpenFOAM Label\"
    FORCE
    )

set(CMAKE_BUILD_TYPE
    \"$ENV{WM_COMPILE_OPTION}\"
    CACHE STRING
    \"Build type (usually Opt, Prof, or Debug)\"
    FORCE
    )
")



file(APPEND ${output_file}
"


# ============================================================================ #
# ---------------------------  Advanced Settings  ---------------------------- #
# ============================================================================ #

include(\${HELYX_PROJECT_DIR}/etc/cmake/defaultSettings/advancedSettings.cmake)

# If you want to overwrite variables set in \"advancedSettings.cmake\", do it
# here.  Don't forget to use the \"FORCE\" keyword.



# ============================================================================ #
# --------------------------- ThirdParty Settings ---------------------------- #
# ============================================================================ #

##include(${HELYX_PROJECT_DIR}/etc/cmake/defaultSettings/thirdPartySettings.cmake)

# NOTE: To consider changing the ThirdParty Settings using the include files
# above. See examples in cmake/cmakeSettingsFiles/exampleLinuxSettings.cmake


# For historical reasons, the following paths are used as default values for
# ThirdParty libraries.  These variables will not appear in the cache.
set(THIRDPARTY_PLATFORMS_PATH_LONG \"\${HELYX_THIRDPARTY_DIR}/platforms/\")
string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG \"\${HELYX_BUILD_PLATFORM}\")
string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG \"\${HELYX_COMPILER_LIB_ARCH}\")
string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG \"\${HELYX_COMPILER_NAME}\")

set(THIRDPARTY_PLATFORMS_PATH_SHORT \"\${THIRDPARTY_PLATFORMS_PATH_LONG}\")
string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG \"\${HELYX_PRECISION_OPTION}\")
string(APPEND THIRDPARTY_PLATFORMS_PATH_LONG \"Int\${HELYX_LABEL_SIZE}\")

")

function(pad_string_with_char string_to_pad character max_length)

    # If centred text, stop one short rather than overunning
    if("${ARGN}" STREQUAL "centre")
        math(EXPR max_length ${max_length}-1)
    endif()

    set(string_to_print "${string_to_pad}")
    string(LENGTH ${string_to_print} string_to_print_length)
    while(${string_to_print_length} LESS ${max_length})
        if("${ARGN}" STREQUAL "centre")
            set(string_to_print "${character}${string_to_print}")
        endif()
        string(APPEND string_to_print "${character}")
        string(LENGTH ${string_to_print} string_to_print_length)
    endwhile()

    if("${ARGN}" STREQUAL "centre")
        math(EXPR max_length ${max_length}+1)
    endif()

    # This is required for centred text, and harmless otherwise
    if(${string_to_print_length} LESS ${max_length})
        string(APPEND string_to_print "${character}")
    endif()

    set(string_to_print "${string_to_print}" PARENT_SCOPE)
endfunction()


set(ENV{PTSCOTCH_ARCH_PATH} "$ENV{SCOTCH_ARCH_PATH}")
set(ENV{PARHIP_ARCH_PATH} "$ENV{KAHIP_ARCH_PATH}")
foreach(library_name
"MPI"
"SCOTCH"
"PTSCOTCH"
"KAHIP"
"PARHIP"
"FFTW"
#
"BOOST"
"TBB"
"CBLOSC"
"OPENVDB"
"OPENCASCADE"
"IMATH"
"ALEMBIC"
"PRECICE"
)
    pad_string_with_char(" ${library_name} " "-" 76 centre)
    set(subtitle "# ${string_to_print} #")
    file(APPEND ${output_file} "\n${subtitle}\n")

    set(compile_flag "$ENV{COMPILE_${library_name}}")
    if("MPI" STREQUAL "${library_name}")
        set(compile_string
            "# MPI is always compiled\n"
        )
    elseif("" STREQUAL "${compile_flag}")
        set(compile_string
            "# ${library_name} will compile if found (${library_name}_REQUIRED not set, defaulting to \"AUTO\")\n"
        )
    else()
        set(compile_string
            "set(${library_name}_REQUIRED\n"
            "    \"${compile_flag}\"\n"
            "    CACHE STRING\n"
            "    \"Trinary flag for forcing compilation of ${library_name} (One of ON, OFF, AUTO)\"\n"
            "    FORCE\n"
            "    )\n"
        )
    endif()
    file(APPEND ${output_file} ${compile_string})

    set(archpath_flag "$ENV{${library_name}_ARCH_PATH}")
    if("" STREQUAL "${archpath_flag}")
        set(archpath_string
            "# ${library_name}_ARCH_PATH not set, system defaults will be used if available\n"
        )
        file(APPEND ${output_file} "${archpath_string}")
    else()
        get_filename_component(real_path "$ENV{${library_name}_ARCH_PATH}" ABSOLUTE)
        normalise_path("${real_path}"
            ${library_name}_ARCH_PATH_STRING
        )
        file(APPEND ${output_file}
            "set(${library_name}_ARCH_PATH\n"
            "    \"${${library_name}_ARCH_PATH_STRING}\"\n"
            "    CACHE PATH\n"
            "    \"Default path on which to look for the ${library_name} library\"\n"
            "    FORCE\n"
            "    )\n"
        )
    endif()

    if("MPI" STREQUAL "${library_name}")
        get_filename_component(FRIENDLY_MPI_NAME "$ENV{MPI_ARCH_PATH}" NAME)
        if ("${FRIENDLY_MPI_NAME}" STREQUAL "$ENV{FOAM_MPI}")
            set(HELYX_MPI_NAME_STRING
                "get_filename_component(FRIENDLY_MPI_NAME \"\${MPI_ARCH_PATH}\" NAME)\n"
                "set(HELYX_MPI_NAME\n"
                "    \"\${FRIENDLY_MPI_NAME}\"\n"
                "    CACHE STRING\n"
                "    \"Human-readable name for MPI, used to set paths for MPI-dependent libraries like Pstream\"\n"
                "    FORCE\n"
                "    )\n"
            )
        else ()
            set(HELYX_MPI_NAME_STRING
                "set(HELYX_MPI_NAME\n"
                "    \"$ENV{FOAM_MPI}\"\n"
                "    CACHE STRING\n"
                "    \"Human-readable name for MPI, used to set paths for MPI-dependent libraries like Pstream\"\n"
                "    FORCE\n"
                "    )\n"
            )
        endif()
        file(APPEND ${output_file} ${HELYX_MPI_NAME_STRING})
    endif()

    set(suffix_flag "$ENV{${library_name}_PATH_SUFFIXES}")
    if("" STREQUAL "${suffix_flag}")
        set(suffix_string
            "# ${library_name}_PATH_SUFFIXES not set.  Only the standard subdirectories (lib,\n"
            "# lib64, inc, include, etc...) of ${library_name}_ARCH_PATH will be searched\n"
        )
    else()
        set(suffix_string
            "set(${library_name}_PATH_SUFFIXES\n"
            "    \"${suffix_flag}\"\n"
            "    CACHE STRING\n"
            "    \"Additional subdirectories of ${archpath_flag} to search\"\n"
            "    FORCE\n"
            "    )\n"
        )
    endif()
    file(APPEND ${output_file} ${suffix_string})

endforeach()


#  As VTK is now detected by RTPP, extra checks are needed
# if VTK_RENDERING_BACKEND or VTK_RENDERING_BACKEND evironment variables are
# explicitly set, it is probably because RTPP is required.
if("$ENV{runTimePostProcessing_REQUIRED}" STREQUAL "")
    set(ENV{runTimePostProcessing_REQUIRED} OFF)
endif()
if(NOT "$ENV{VTK_ARCH_PATH}")
    set(ENV{VTK_ARCH_PATH} "\${HELYX_PROJECT_DIR}/../../EVTK/ext")
else()
    set(ENV{runTimePostProcessing_REQUIRED} AUTO)
endif()
if(NOT "$ENV{VTK_RENDERING_BACKEND}")
    set(ENV{VTK_RENDERING_BACKEND} OSMesa)
    else()
    set(ENV{runTimePostProcessing_REQUIRED} AUTO)
endif()

pad_string_with_char(" RTPP - EVTK " "-" 76 centre)
set(subtitle "# ${string_to_print} #")
file(APPEND ${output_file} "\n${subtitle}\n")

file(APPEND ${output_file}
"
# ============================================================================ #
# ---------------------------- EVTK/RTPP Settings ---------------------------- #
# ============================================================================ #

# Engys fork of VTK (EVTK) is a dependency of Run-Time Post-Processing (RTPP).
# Therefore, RTPP will only compile only if EVTK is found.
# By default, RTPP compilation is disabled. To enable RTPP compilation set
# runTimePostProcessing_REQUIRED to \"ON\" or \"AUTO\". When enabled, RTPP will look
# for EVTK and if the requirements are met, RTPP will compile.
#
# runTimePostProcessing_REQUIRED not set: defaults to \"AUTO\";
#
# VTK_REQUIRED not set: defaults to \"AUTO\"; it will search for EVTK using the
# VTK_ARCH_PATH, and the HELYX GUI EVTK path ( ~/Engys/HELYX/**/GUI/ext/ );
#
# VTK_ARCH_PATH not set: only HELYX GUI EVTK path ( ~/Engys/HELYX/**/GUI/ext/ )
# will be searched and used if available
#
# VTK_PATH_SUFFIXES not set:  Only the standard subdirectories
# (lib, lib64, inc, include, etc...) of VTK_ARCH_PATH will be searched.
#
# For custom EVTK compilation, uncomment the lines below, and set the
# VTK_ARCH_PATH to the compiled EVTK

set(runTimePostProcessing_REQUIRED
    $ENV{runTimePostProcessing_REQUIRED}
    CACHE STRING
    \"To compile RTPP. It requires EVTK\"
    FORCE
    )

set(VTK_ARCH_PATH
    $ENV{VTK_ARCH_PATH}
    CACHE PATH
    \"Default path on which to look for the VTK library\"
    FORCE
    )
set(VTK_RENDERING_BACKEND
    $ENV{VTK_RENDERING_BACKEND}
    CACHE STRING
    \"The VTK rendering backend (one of OSMesa, OpenGL2)\"
    FORCE
    )
")

file(APPEND ${output_file} "\n")


file(APPEND ${output_file}
"
# ============================================================================ #
# Set sensible defaults for ThirdParty compile switches, etc...
include(\${HELYX_PROJECT_DIR}/etc/cmake/parseThirdPartySettings.cmake)

")

get_filename_component(output_file_abs "${output_file}" ABSOLUTE)
string(CONCAT s
    "userSettings.cmake was automatically generated from the current "
    "environment variables.\n\n"
    "\tThe CMake build system does not use environment variables to configure"
    " HELYX.\n"
    "\tAll future changes to the HELYX configuration should be made in the "
    "following file:\n"
    "\t\t${output_file_abs}\n\n"
)
message(STATUS "${s}")
