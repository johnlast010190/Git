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
    If hpc module variable HELYX_HPC_MODULE is defined in userSettingsFile.cmake
    then, a module file for Helyx is generated based on current cmake configuration

[----------------------------------------------------------------------------]]

# Re-loading HELYX_RUNTIME_OUTPUT_DIRECTORY and HELYX_LIBRARY_OUTPUT_DIRECTORY
if("${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/bin")
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/bin")
else()
    set(HELYX_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endif()

if("${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" STREQUAL "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib")
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/lib")
else()
    set(HELYX_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
endif()

# ------------------------------ Setting Strings ----------------------------- #

set(HPC_MODULE_HEADER
"#%Module1.0
#
${CONFIG_FILE_HEADER}
#
# Description
#     HELYX environment module file for HPCs.
#     This is an automatically-generated environment module of HELYX.
#     This file exports CMake variables as environment variables, which define
#     a viable environment to run HELYX.
#
#     This environment-module of HELYXcore is on experimental state.

#
# Comments
#     Copy the content of 'environmentModule/' to the system moduleFile directory,
#     i.e., \"HELYXcore/${HELYX_PROJECT_VERSION}-${HELYX_OPTIONS}\"
#
#     Variables of this moduleFile are defined based on their absolute path.
#     If binaries were provided instead of building HELYX for the current system,
#     please use this module file as template and edit the variables manually.
#
#     If dependency modules are required, they should be manually included with:
#        module load modulefile modulefile
#
#     If system MPI is being used, remove MPI paths from $PATH and $LD_LIBRARY_PATH
#
#
# Exemples of 'environment-modules' initialization
#     C Shell initialization (and derivatives):
#        source /usr/share/Modules/init/csh
#        module load modulefile modulefile ...
#
#     Bourne Shell (sh) (and derivatives):
#        . /usr/share/Modules/init/sh
#        module load modulefile modulefile ...
#
#     Python:
#        import os
#        exec(open('/usr/share/Modules/init/python.py').read())
#        module('load', 'modulefile', 'modulefile', '...')
#
# [----------------------------------------------------------------------------]

proc ModulesHelp{ } {
    global-version HELYXcore-${HELYX_PROJECT_VERSION}
    puts stderr \"\\thelyx - loads the HELYX cmake environment\"
    puts stderr \"\\n\\tVersion $version\\n\"
}

module-whatis \"Loads the HELYXcore-${HELYX_PROJECT_VERSION} environment\"

if [ module-info mode load ] {
    puts stderr \"
Module 'HELYXcore/${HELYX_PROJECT_VERSION}-${HELYX_OPTIONS}' loaded.

    This environment-module of HELYXcore is on experimental state
    and variables are defined based on their absolute path.

    If binaries were provided instead of building HELYX for the current system,
    please use this module file as a template and edit the variables manually.

    If system MPI is being used, make sure MPI paths are not defined in
    PATH and LD_LIBRARY_PATH variables of the module file.\n\"
}

")

string(REPLACE "#!/bin/bash"
    "#"
    HPC_MODULE_HEADER  # Output var
    "${HPC_MODULE_HEADER}"
    )

set(HPC_MODULE_PREPEND_PATH
"# Prepend paths

# Explicitly defining HELYX_THIRDPARTY_DIR for MPI related variables
set HELYX_THIRDPARTY_DIR \"${HELYX_PROJECT_DIR}/../ThirdParty-${HELYX_THIRDPARTY_VERSION}\"

prepend-path PATH \"${HELYX_PROJECT_DIR}/bin\"
prepend-path PATH \"${HELYX_RUNTIME_OUTPUT_DIRECTORY}\"
prepend-path PATH \"${HELYX_PROJECT_DIR}/wmake\"
prepend-path PATH \"${MPIEXEC_PATH}\"

prepend-path LD_LIBRARY_PATH \"${MPI_LIBRARY_PATHS_STRING}\"
")

set(HPC_MODULE_EXTRA_MODULES "# Dependency modules\n")

set(HPC_MODULE_SETENV
"# Environmental variables

setenv HELYX_SIGFPE \"${HELYX_SIGFPE}\"

setenv HELYX_PROJECT_DIR \"${HELYX_PROJECT_DIR}\"

setenv HELYX_SETTINGS_FILE \"${HELYX_SETTINGS_FILE}\"
setenv HELYX_OPTIONS \"${HELYX_OPTIONS}\"

setenv HELYX_SRC \"${HELYX_SRC}\"
setenv HELYX_ETC \"${HELYX_PROJECT_DIR}/etc\"
setenv HELYX_MODULES \"${HELYX_MODULES}\"
setenv HELYX_APPLICATIONS \"${HELYX_APPLICATIONS}\"
setenv HELYX_SOLVERS \"${HELYX_SOLVERS}\"
setenv HELYX_UTILITIES \"${HELYX_UTILITIES}\"
setenv HELYX_TUTORIALS \"${HELYX_PROJECT_DIR}/examples\"

setenv HELYX_RUNTIME_OUTPUT_DIRECTORY \"${HELYX_RUNTIME_OUTPUT_DIRECTORY}\"
setenv HELYX_LIBRARY_OUTPUT_DIRECTORY \"${HELYX_LIBRARY_OUTPUT_DIRECTORY}\"

setenv HELYX_CONFIG \"${HELYX_PROJECT_DIR}/etc/dictData\"

setenv HELYX_THIRDPARTY_VERSION \"${HELYX_THIRDPARTY_VERSION}\"
# setenv HELYX_THIRDPARTY_DIR \"${HELYX_THIRDPARTY_DIR}\"
setenv HELYX_THIRDPARTY_DIR \"${HELYX_PROJECT_DIR}/../ThirdParty-${HELYX_THIRDPARTY_VERSION}\"

setenv HELYX_MPI_NAME \"${HELYX_MPI_NAME}\"
setenv MPI_ARCH_PATH \"${MPI_DIR_STRING}\"
setenv OPAL_PREFIX \"${MPI_DIR_STRING}\"

setenv HELYX_BUILD_PLATFORM \"${HELYX_BUILD_PLATFORM}\"
setenv HELYX_COMPILER_LIB_ARCH \"${HELYX_COMPILER_LIB_ARCH}\"
setenv HELYX_COMPILER_NAME \"${HELYX_COMPILER_NAME}\"
setenv HELYX_PRECISION_OPTION \"${HELYX_PRECISION_OPTION}\"
setenv HELYX_LABEL_SIZE \"${HELYX_LABEL_SIZE}\"
setenv HELYX_BUILD_TYPE \"${HELYX_BUILD_TYPE}\"
setenv HELYX_PROJECT_VERSION \"${HELYX_PROJECT_VERSION}\"
")

## An alternative is to source 'active' build directly instead of exporting the variables
#set(HPC_MODULE_FOOTER
#"
## HELYX sourceable files
#switch -- [module-info shelltype]{
#    sh{
#        source-sh bash $HELYX_PROJECT_DIR/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc
#    }
#    csh{
#        source-sh csh $HELYX_PROJECT_DIR/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.cshrc
#    }
#}
#")

set(output_string "")
string(CONCAT output_string "${output_string}"
"${HPC_MODULE_HEADER}

# [----------------------------------------------------------------------------]

${HPC_MODULE_EXTRA_MODULES}

# [----------------------------------------------------------------------------]

${HPC_MODULE_PREPEND_PATH}

# [----------------------------------------------------------------------------]

${HPC_MODULE_SETENV}

# [----------------------------------------------------------------------------]
")

# Whatever else ends up in the config file, it must end in a newline to be POSIX
string(CONCAT output_string "${output_string}" "\n")

# -------------------------- Writting HPC module file ------------------------ #

set(HPC_MODULE_FILE ${CMAKE_BINARY_DIR}/extraSourceableFiles/environmentModule/HELYXcore/${HELYX_PROJECT_VERSION}-${HELYX_OPTIONS})

file(WRITE ${HPC_MODULE_FILE} ${output_string})

# message(STATUS "Generated a HPC environment-module file at the following location:
#     ${HPC_MODULE_FILE}")
