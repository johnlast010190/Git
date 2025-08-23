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
    This file contains strings that are common to two or more of the sourceable
    file generators.

[----------------------------------------------------------------------------]]


set(CONFIG_FILE_HEADER
"#!/bin/bash
# [----------------------------------------------------------------------------]
# |       o        |
# |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
# |   o   O   o    |  Version : ${HELYX_PROJECT_VERSION}
# |    o     o     |  ENGYS Ltd. <http://engys.com/>
# |       o        |
# [----------------------------------------------------------------------------]
# License
#     This file is part of HELYXcore.
#     HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
#
#     HELYXcore is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     HELYXcore is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright
#     (c) 2020-2021 Engys Ltd"
)


get_filename_component(REAL_TP_DIR "${HELYX_THIRDPARTY_DIR}" REALPATH)
file(RELATIVE_PATH RELATIVE_TP_DIR ${HELYX_PROJECT_DIR} ${HELYX_THIRDPARTY_DIR})
string(REPLACE
    "ThirdParty-${HELYX_THIRDPARTY_VERSION}"
    ThirdParty-\${HELYX_THIRDPARTY_VERSION}
    TP_DIR_STRING
    "${RELATIVE_TP_DIR}"
    )
set(TP_DIR_STRING "\${HELYX_PROJECT_DIR}/${TP_DIR_STRING}")


if("${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/bin")
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "\${HELYX_PROJECT_DIR}/platforms/\${HELYX_OPTIONS}/bin")
else()
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endif()

if("${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib")
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "\${HELYX_PROJECT_DIR}/platforms/\${HELYX_OPTIONS}/lib")
else()
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
endif()

get_filename_component(REAL_HELYX_PROJECT_DIR "${HELYX_PROJECT_DIR}" REALPATH)
string(REPLACE
    "${REAL_HELYX_PROJECT_DIR}" "\${HELYX_PROJECT_DIR}"
    HELYX_SETTINGS_FILE_STRING
    "${HELYX_SETTINGS_FILE}"
    )

set(PATH_VARIABLES_STRING
"# -----------------------------  Path variables  ----------------------------- #

# It's impossible to find the path of a sourced script in POSIX, so we must use
# specialised methods from the derived shells.  If the provided bash and zsh
# methods fail, then HELYX_PROJECT_DIR falls back to the absolute location at
# compile-time
script_dir=\"\${BASH_SOURCE:-\${ZSH_NAME:+$0}}\"
[ -n \"$script_dir\" ] && \\
HELYX_PROJECT_DIR=$(\\cd $(dirname $script_dir)/.. && \\pwd -L) || \\
HELYX_PROJECT_DIR=\"${HELYX_PROJECT_DIR}\"
export HELYX_PROJECT_DIR
unset script_dir

export HELYX_SETTINGS_FILE=\"${HELYX_SETTINGS_FILE_STRING}\"

export HELYX_OPTIONS=\"${HELYX_OPTIONS}\"

export HELYX_SRC=\"\${HELYX_PROJECT_DIR}/src\"
export HELYX_ETC=\"\${HELYX_PROJECT_DIR}/etc\"
export HELYX_MODULES=\"\${HELYX_PROJECT_DIR}/modules\"
export HELYX_APPLICATIONS=\"\${HELYX_PROJECT_DIR}/applications\"
export HELYX_SOLVERS=\"\${HELYX_APPLICATIONS}/solvers\"
export HELYX_UTILITIES=\"\${HELYX_APPLICATIONS}/utilities\"
export HELYX_TUTORIALS=\"\${HELYX_PROJECT_DIR}/examples\"

export HELYX_RUNTIME_OUTPUT_DIRECTORY=\"${HELYX_RUNTIME_OUTPUT_DIRECTORY}\"
export HELYX_LIBRARY_OUTPUT_DIRECTORY=\"${HELYX_LIBRARY_OUTPUT_DIRECTORY}\"

export HELYX_CONFIG=\"\${HELYX_PROJECT_DIR}/etc/dictData\"
export FOAM_CONFIG=\"\${HELYX_PROJECT_DIR}/etc/dictData\"

export HELYX_THIRDPARTY_VERSION=\"${HELYX_THIRDPARTY_VERSION}\"
export HELYX_THIRDPARTY_DIR=${TP_DIR_STRING}")


set(HELPER_VARIABLES_STRING
"# ------------------------------  IDE variables  ----------------------------- #
# The following variables are useful when configuring IDE's.  They are NOT used
# in the compilation steps.  These settings must be changed in
# userSettings.cmake (or equivalent).

export HELYX_BUILD_PLATFORM=${HELYX_BUILD_PLATFORM}
export HELYX_COMPILER_LIB_ARCH=${HELYX_COMPILER_LIB_ARCH}
export HELYX_COMPILER_NAME=${HELYX_COMPILER_NAME}
export HELYX_PRECISION_OPTION=${HELYX_PRECISION_OPTION}
export HELYX_LABEL_SIZE=${HELYX_LABEL_SIZE}
export HELYX_BUILD_TYPE=${CMAKE_BUILD_TYPE}
export HELYX_PROJECT_VERSION=${HELYX_PROJECT_VERSION}

# Useful variables for compilation of ThirdParty libraries
# (to make sure we compile ThirdParty with same compiler used for Helyx)
export HELYX_CC=${CMAKE_C_COMPILER}
export HELYX_CXX=${CMAKE_CXX_COMPILER}
export THIRDPARTY_INSTALL_DIR=${TP_DIR_STRING}/platforms/${HELYX_BUILD_PLATFORM}${HELYX_COMPILER_LIB_ARCH}${HELYX_COMPILER_NAME}


# ---------------------------------  Aliases  -------------------------------- #
alias src='cd \${HELYX_SRC}'
alias run='cd \${HELYX_RUN}'
alias tut='cd \${HELYX_TUTORIALS}'
alias sol='cd \${HELYX_SOLVERS}'
alias app='cd \${HELYX_APPLICATIONS}'
alias util='cd \${HELYX_UTILITIES}'
alias foam='cd \${HELYX_PROJECT_DIR}'
alias helyx='cd \${HELYX_PROJECT_DIR}'


# ---------------------------------  Misc  -------------------------------- #
# source emake auto completion
emakeAutocompleteFile=$HELYX_PROJECT_DIR/platforms/$HELYX_OPTIONS/.emakeAutocomplete.sh
if [[ -f $emakeAutocompleteFile ]] && [[ -x $emakeAutocompleteFile ]]; then
    # for auto-completion to work with zsh
    [ -n \"\$ZSH_VERSION\" ] && autoload -U bashcompinit && bashcompinit
    . $emakeAutocompleteFile
fi

# source auto completion only for bash ge 4.2
bashCompletionFile=$HELYX_PROJECT_DIR/etc/config.sh/bash_completion
[ -n \"\$BASH_VERSION\" ] && bash_version=$(echo \"\$BASH_VERSION\" | perl -pe '(\$_)=/([0-9]+[.][0-9]+)/')
if [[ -f $bashCompletionFile ]] && [[ -n \"\$bash_version\" ]] && \
[[ \"\$(printf '%s\\n' \"4.2\" \"$bash_version\" | sort -V | head -n1)\" == \"4.2\" ]]; then
    . $bashCompletionFile
fi
")

set(HELPER_MPI_VARIABLES_STRING
"# ------------------------------  MPI variables  ----------------------------- #
# MPI variables are needed to be able to make the code properly MPI dependent
# This variables will be also defined later in the complete sourceable file
export HELYX_MPI_NAME=${HELYX_MPI_NAME}
export MPI_ARCH_PATH=${MPI_ARCH_PATH}

")
