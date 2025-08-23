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



# ============================================================================ #
# ---------------------------- activeBuild.cshrc ----------------------------- #
# ============================================================================ #

# Always populate extraSourceableFiles dir to make packing reproducable
set(ACTIVE_BUILD_FILE ${CMAKE_BINARY_DIR}/extraSourceableFiles/activeBuild-${HELYX_PRECISION_OPTION}.cshrc)

set(active_build_output_string
"${CONFIG_FILE_HEADER}
#
# Description
#     This is an automatically-generated file that points to the latest rc file.
#     It's re-generated whenever userSettingsConfig.cmake is run.
#
#     Any changes made to this file will be lost when userSettingsConfig.cmake
#     is re-run, e.g. every time the CMakeCache is refreshed.
#
# [----------------------------------------------------------------------------]

echo Error:  activeBuild is only implemented for POSIX shells.  Please manually source the following file:
echo '    $HELYX_PROJECT_DIR/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.cshrc'
exit 1
")

string(REPLACE "#!/bin/bash"
    "#!/bin/csh"
    active_build_output_string  # Output var
    "${active_build_output_string}"
    )

# Whatever else ends up in the config file, it must end in a newline to be POSIX
# compliant
string(CONCAT active_build_output_string "${active_build_output_string}" "\n")
file(WRITE ${ACTIVE_BUILD_FILE} ${active_build_output_string})



# ============================================================================ #
# -------------------------- ${HELYX_OPTIONS}.cshrc -------------------------- #
# ============================================================================ #

set(SH_FILE ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc)
# Always populate extraSourceableFiles dir to make packing reproducable
set(CSH_FILE ${CMAKE_BINARY_DIR}/extraSourceableFiles/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.cshrc)

# Always generate shrc files, because they're needed by emake
if(EXISTS "${SH_FILE}")
    file(READ "${SH_FILE}" output_string)
else()
    message(WARN
        "POSIX sourceable file not found, unable to generate csh sourceable file"
        )
    set(output_string "")
endif()

# string(REGEX REPLACE <regular_expression>
#        <replacement_expression> <output_variable>
#        <input> [<input>...])

# Replace line continuations in PATH and LD_LIBRARY_PATH
string(REGEX REPLACE ":\\\\\n"
    "\\\\:"
    output_string  # Output var
    "${output_string}"
    )

string(REPLACE "#!/bin/bash"
    "#!/bin/csh"
    output_string  # Output var
    "${output_string}"
    )

# Replace "" PATH and LD_LIBRARY_PATH
string(REGEX REPLACE "export PATH=\"([^\"]*)\""
"if (! $?PATH) then
setenv PATH
endif
setenv PATH \\1"
    output_string  # Output var
    "${output_string}"
    )
string(REGEX REPLACE "export LD_LIBRARY_PATH=\"([^\"]*)\""
"if (! $?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH
endif
setenv LD_LIBRARY_PATH \\1"
    output_string  # Output var
    "${output_string}"
    )

string(REPLACE
# Replace this:
"# It's impossible to find the path of a sourced script in POSIX, so we must use
# specialised methods from the derived shells.  If the provided bash and zsh
# methods fail, then HELYX_PROJECT_DIR falls back to the absolute location at
# compile-time
script_dir=\"\${BASH_SOURCE:-\${ZSH_NAME:+$0}}\"
[ -n \"$script_dir\" ] && \\
HELYX_PROJECT_DIR=$(\\cd $(dirname $script_dir)/.. && \\pwd -L) || \\
HELYX_PROJECT_DIR=\"${HELYX_PROJECT_DIR}\"
export HELYX_PROJECT_DIR
unset script_dir"
# With this:
"setenv HELYX_PROJECT_DIR `lsof +p $$ |& sed -n -e 's@[^/]*@@' -r -n -e 's@(/[^/]*)(/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.cshrc).*@\\1@p'`"
    output_string  # Output var
    "${output_string}"
    )

# Replace export statements
string(REGEX REPLACE "export ([^=]+)=([^\n]*)"
    "setenv \\1 \\2"
    output_string  # Output var
    "${output_string}"
    )

# Replace aliases
string(REGEX REPLACE "alias ([^=]+)=([^\n]*)"
    "alias \\1 \\2"
    output_string  # Output var
    "${output_string}"
    )


# Remove autocomplete
string(REGEX REPLACE
"# ---------------------------------  Misc  -------------------------------- #
.*
.*
.*
.*
"
    ""
    output_string  # Output var
    "${output_string}"
    )

# Whatever else ends up in the config file, it must end in a newline to be POSIX
# compliant
string(CONCAT output_string "${output_string}" "\n")
file(WRITE ${CSH_FILE} ${output_string})

if ("csh" IN_LIST RUNTIME_SHELL)
    file(COPY "${ACTIVE_BUILD_FILE}" "${CSH_FILE}"
        DESTINATION "${HELYX_PROJECT_DIR}/platforms/"
    )
    configure_file(${HELYX_PROJECT_DIR}/etc/cmake/templates/template-unsetActiveBuild.cshrc
        ${HELYX_PROJECT_DIR}/platforms/unsetActiveBuild.cshrc
        @ONLY
    )
endif()
