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

[----------------------------------------------------------------------------]]


set(MPI_IN_THIRDPARTY FALSE)
if(NOT "" STREQUAL "${MPI_ARCH_PATH}")
    get_filename_component(REAL_MPI_DIR "${MPI_ARCH_PATH}" REALPATH)
    string(REPLACE
        ${REAL_TP_DIR} \${HELYX_THIRDPARTY_DIR}
        MPI_DIR_STRING ${REAL_MPI_DIR}
        )
    if(NOT ${MPI_DIR_STRING} STREQUAL ${REAL_MPI_DIR})
        set(MPI_IN_THIRDPARTY TRUE)
    endif()
else()
    set(REAL_MPI_DIR "")
    set(MPI_DIR_STRING "\"\"")
endif()

set(output_string "")
string(CONCAT output_string "${output_string}"
"${CONFIG_FILE_HEADER}
#
# Description
#     This is an automatically-generated file that exports CMake variables as
#     environment variables.  It also appends some paths to the LD_LIBRARY_PATH.
#     It's re-generated whenever userSettingsConfig.cmake is run.
#
#     These variables are exported to the shell for convenience only, and have
#     no impact on the compilation environment.
#
#     Any changes made to this file will be lost when userSettingsConfig.cmake
#     is re-run, e.g. every time the CMakeCache is refreshed.
#
# [----------------------------------------------------------------------------]


# ---------------------------  Run-time variables  --------------------------- #

#  Switch to catch floating-point errors (set by default)
export HELYX_SIGFPE=

#  Switch to initialise memory to NaN (unset by default)
#export HELYX_SETNAN=


${PATH_VARIABLES_STRING}

")

string(CONCAT output_string "${output_string}"
"
export HELYX_MPI_NAME=${HELYX_MPI_NAME}  # Required for dynamicCode Makefile template
export MPI_ARCH_PATH=${MPI_DIR_STRING}
")

# Set OPAL_PREFIX only if MPI is in ThirdParty. This makes MPI from ThirdParty
# re-locatable.  Anything else risks being unexpected, and breaking users
# existing external MPI's.
if("${MPI_IN_THIRDPARTY}" STREQUAL TRUE)
    string(CONCAT output_string "${output_string}"
    "export OPAL_PREFIX=${MPI_DIR_STRING}\n")
endif()

# Add the path to mpiexec to PATH
if(NOT MPI_AWARE_COMPILER)

    if(EXISTS ${MPIEXEC})
        get_filename_component(MPIEXEC_PATH "${MPIEXEC}" DIRECTORY)
    elseif(NOT "" STREQUAL "${MPI_ARCH_PATH}" AND NOT "${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
        message(WARNING
            "Could not determine path to mpiexec for export, using MPI_ARCH_PATH/bin instead.
Please check that MPI is configured correctly at run-time.
        ")
        set(MPIEXEC_PATH "${MPI_ARCH_PATH}/bin")
    elseif(NOT "" STREQUAL "${MPI_ARCH_PATH}" AND "${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
        message(WARNING
            "Its not possible to determine path to mpiexec for Windows cross-compilation exports.
mpiexec is not installed in Windows ThirdParty.
Using hard coded path \"C:\\\\Program Files\\\\Microsoft MPI\\\\Bin\" in the *.bat sourceable file.
    ")
    set(MPIEXEC_PATH "")
    else()
        message(WARNING "MPI found, but MPI_ARCH_PATH set to \"${MPI_ARCH_PATH}\"
Cannot prepend path to mpiexec to PATH.  It is up to you to make sure your \
mpi environment is configured correctly at run-time.
        ")
        set(MPIEXEC_PATH "")
    endif()

    get_filename_component(MPIEXEC_PATH "${MPIEXEC_PATH}" REALPATH)
    string(REPLACE
        "${REAL_TP_DIR}" \${HELYX_THIRDPARTY_DIR}
        MPIEXEC_PATH "${MPIEXEC_PATH}"
        )

    string(CONCAT output_string "${output_string}"
    "
# Now add directories to PATH, including path to the mpiexec for the MPI
# against which HELYX is compiled:"
    )

else()

    string(CONCAT output_string "${output_string}"
    "
# Now add directories to PATH. MPI not added to path since this is a
# MPI-aware compile which is already exposed in the system:"
    )

endif()

# Exporting to PATH
string(CONCAT output_string "${output_string}"
"
export PATH=\"\\
\${HELYX_PROJECT_DIR}/bin:\\
\${HELYX_RUNTIME_OUTPUT_DIRECTORY}:\\
\${HELYX_PROJECT_DIR}/wmake:\\
")

# Don't mess with system paths!
# Would be elegant to do this with CMake variables, but CMake doesn't appear to
# expose this information...
# Note that this "bug" has only ever been seen when cross-compiling for Windows,
# but messing with where system paths appears on PATH seems like a bad idea
# anyway.
if( NOT "" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "MSwindows" STREQUAL "${HELYX_SYSTEM_NAME}" AND
    NOT "/bin" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "/usr/bin" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "/usr/local/bin" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "/sbin" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "/usr/sbin" STREQUAL "${MPIEXEC_PATH}" AND
    NOT "/usr/local/sbin" STREQUAL "${MPIEXEC_PATH}"
)
    string(CONCAT output_string "${output_string}"
        "${MPIEXEC_PATH}:\\
")
endif()

if(NOT "" STREQUAL "${CMAKE_TEST_OUTPUT_DIRECTORY}")
    string(CONCAT output_string "${output_string}"
        "\${HELYX_PROJECT_DIR}/platforms/\${HELYX_OPTIONS}/test:\\
")
endif()

string(CONCAT output_string "${output_string}"
"$PATH\"

")

if(NOT MPI_AWARE_COMPILER)

    foreach(lib IN LISTS MPI_LIBRARY_PATHS)
        get_filename_component(lib "${lib}" REALPATH)
        list(APPEND RESOLVED_MPI_LIBRARY_PATHS "${lib}")
    endforeach()

    # Add thirdParty libraries that depend on MPI to LD_LIBRARY_PATH,
    # but only if the packages from ThirdParty are used
    foreach(lib IN LISTS PARALLEL_THIRDPARTY_LIBRARIES)
        set(lib_path "${${lib}_ARCH_PATH}")
        get_filename_component(lib_path "${lib_path}" REALPATH)
        if(${lib}_FOUND AND "${lib_path}" MATCHES "${HELYX_THIRDPARTY_DIR}" )
            list(APPEND RESOLVED_MPI_LIBRARY_PATHS "${lib_path}/lib/${HELYX_MPI_NAME}")
        endif()
    endforeach()

    string(REPLACE
        ";" ":\\\n"
        MPI_LIBRARY_PATHS_STRING "${RESOLVED_MPI_LIBRARY_PATHS}"
        )
    string(REPLACE
        "${REAL_TP_DIR}" \${HELYX_THIRDPARTY_DIR}
        MPI_LIBRARY_PATHS_STRING "${MPI_LIBRARY_PATHS_STRING}"
        )

    string(CONCAT output_string "${output_string}"
    "# The MPI libraries must be appended to LD_LIBRARY_PATH for MPI to be re-
# locatable.  Again, use the hard-coded compile-time value and include the
# default thirdParty path afterwards."
    )

else()

    string(CONCAT output_string "${output_string}"
    "# No need to appended MPI libraries to LD_LIBRARY_PATH since this is a
# MPI-aware compile which is already exposed by the system:"
    )

endif()

string(CONCAT output_string "${output_string}"
"# No need to add HELYX libraries, as we now use RUNPATH.
export LD_LIBRARY_PATH=\"\\
${MPI_LIBRARY_PATHS_STRING}:\\
")

string(CONCAT  output_string "${output_string}"
"\${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib/${HELYX_MPI_NAME}:\\
")

string(CONCAT output_string "${output_string}"
"$LD_LIBRARY_PATH\"
")

# This is a deliberately undocumented feature, implemented as a work-around for
# adding extra dependencies to extremely old Linux platforms.  Don't use this
# unless you have to - if users have system-specific requirements, they're
# generally expected to look after them themselves.
if(NOT "" STREQUAL "${HELYX_ADDITIONAL_LD_LIBRARY_PATHS}")
    string(CONCAT output_string "${output_string}"
        "\nexport LD_LIBRARY_PATH=\"${HELYX_ADDITIONAL_LD_LIBRARY_PATHS}:$LD_LIBRARY_PATH\"\n"
    )
endif()

string(CONCAT
    output_string "${output_string}"
    "\n\n${HELPER_VARIABLES_STRING}"
    )


# Warn if ThirdParty GCC/G++ is being used
set(third_party_gcc_warning
"
# ---------------------------------  Warning  -------------------------------- #
echo -e \"\nWARNING:  Using custom ThirdParty Gcc and G++.
Please make sure that relevant compiler libraries and \
dependencies (gmp, mpc, mpfr) are listed in LD_LIBRARY_PATH. \n\"
")

if("${third_party_gcc}" STREQUAL "TRUE")
    string(CONCAT
        output_string "${output_string}"
        "\n\n${third_party_gcc_warning}"
        )
endif()

# Warn if old environment variables (WM_* and FOAM_*) are still defined
# This check must be done explicitly here - the environment can change without
# triggering a re-compilation of this file!
# message(FATAL_ERROR "ENVIRONMENT_VARIABLE_WARNING_WHITELIST = ${ENVIRONMENT_VARIABLE_WARNING_WHITELIST}")

set(whitelist_string "whitelist=''")
if(NOT "" STREQUAL "${ENVIRONMENT_VARIABLE_WARNING_WHITELIST}")
    list(GET ENVIRONMENT_VARIABLE_WARNING_WHITELIST 0 search_string)
    set(whitelist_string "whitelist='${search_string}'")
    list(REMOVE_AT ENVIRONMENT_VARIABLE_WARNING_WHITELIST 0)
    foreach(search_string ${ENVIRONMENT_VARIABLE_WARNING_WHITELIST})
        string(CONCAT whitelist_string
            "${whitelist_string}\n"
            "whitelist+='|${search_string}'"
        )
    endforeach()
endif()

string(CONCAT
output_string "${output_string}"
"


# ---------------------------  Environment check  ---------------------------- #

# These variables from the ENVIRONMENT_VARIABLE_WARNING_WHITELIST CMake cache
# variable.  For more details, see
# \"\${HELYX_PROJECT_DIR}/etc/cmake/defaultSettings/advancedSettings.cmake\".
${whitelist_string}

declare -a suspicious_variables=()
for env_var in $(env); do
    [[ \"$env_var\" =~ '^FOAM_|^WM_'  ]] && \\
    [[ ! $whitelist || ! \"$env_var\" =~ \${whitelist} ]] && \\
    suspicious_variables+=(\"$env_var\")
done

if [ -n \"$suspicious_variables\" ]; then
    printf '%s\\n' \\
    'WARNING:  Some environment variables prefixed with \"FOAM_\" or \"WM_\" '\\
    'were detected.  These variables will have no effect on the HELYXcore '\\
    'compilation but may have undesirable effects on the HELYX run-time '\\
    'environment.'\\
    'The variables detected are as follows:'
    printf '\\t%s\\n' \"\${suspicious_variables[@]}\"
    printf '%s\\n' \\
    '\n(In order to silence this warning, either unset these variables or '\\
    'add them to the \"ENVIRONMENT_VARIABLE_WARNING_WHITELIST\" in your user '\\
    'settings file).'
fi


# ---------------------------------------------------------------------------- #
")

# Whatever else ends up in the config file, it must end in a newline to be POSIX
# compliant
string(CONCAT output_string "${output_string}" "\n")

file(WRITE "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc" "${output_string}")

configure_file(${HELYX_PROJECT_DIR}/etc/cmake/templates/template-unsetActiveBuild.shrc
    ${HELYX_PROJECT_DIR}/platforms/unsetActiveBuild.shrc
    @ONLY
)