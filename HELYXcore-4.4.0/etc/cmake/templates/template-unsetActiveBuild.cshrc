# [---------------------------------------------------------------------------]
# |       o        |
# |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
# |   o   O   o    |  Version : 4.4.0
# |    o     o     |  ENGYS Ltd. <http://engys.com/>
# |       o        |
# [---------------------------------------------------------------------------]
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

# Copyright
#     (C) 2011-2016 OpenFOAM Foundation
#     (C) 2016-2022 OpenCFD Ltd.
#     (C) 2023 Engys Ltd.
#
# File
#     platforms/unsetActiveBuild.chrc
#     - need to be sourced manually
#
# Description
#     Clear as many HELYX environment settings as possible
#
#------------------------------------------------------------------------------

# Clean standard environment variables (PATH, LD_LIBRARY_PATH, MANPATH)

unset foamClean
if ( $?HELYX_PROJECT_DIR ) then
    set foamClean="$HELYX_PROJECT_DIR/bin/foamCleanPath"
    if ( ! -f "$foamClean" || ! -x "$foamClean" ) then
    echo "WARNING: Not cleaning PATH and LD_LIBRAY_PATH as \
'$HELYX_PROJECT_DIR/bin/foamCleanPath' requires execute permissions"
    unset foamClean
    endif
endif

# The old dirs to be cleaned from the environment variables
#set foamOldDirs="$HELYX_PROJECT_DIR $HELYX_THIRD_PARTY_DIR $HOME/$HELYX_PROJECT/$LOGNAME $FOAM_SITE_APPBIN $FOAM_SITE_LIBBIN"
set foamOldDirs="${HELYX_PROJECT_DIR} ${HELYX_PROJECT_DIR}/bin \
${HELYX_RUNTIME_OUTPUT_DIRECTORY} ${HELYX_LIBRARY_OUTPUT_DIRECTORY} \
${HELYX_PROJECT_DIR}/wmake ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/test \
@MPIEXEC_PATH@ @MPI_LIBRARY_PATHS_STRING@ \
@HELYX_ADDITIONAL_LD_LIBRARY_PATHS@"



#------------------------------------------------------------------------------
# Unset HELYX_* environment variables

# Path variables
unsetenv HELYX_PROJECT_DIR
unsetenv HELYX_SETTINGS_FILE
unsetenv HELYX_OPTIONS
unsetenv HELYX_SRC
unsetenv HELYX_ETC
unsetenv HELYX_MODULES
unsetenv HELYX_APPLICATIONS
unsetenv HELYX_SOLVERS
unsetenv HELYX_UTILITIES
unsetenv HELYX_TUTORIALS
unsetenv HELYX_RUNTIME_OUTPUT_DIRECTORY
unsetenv HELYX_LIBRARY_OUTPUT_DIRECTORY
unsetenv HELYX_CONFIG
unsetenv HELYX_THIRDPARTY_VERSION
unsetenv HELYX_THIRDPARTY_DIR
unsetenv HELYX_SIGFPE

# IDE variables
unsetenv HELYX_BUILD_PLATFORM
unsetenv HELYX_COMPILER_LIB_ARCH
unsetenv HELYX_COMPILER_NAME
unsetenv HELYX_PRECISION_OPTION
unsetenv HELYX_LABEL_SIZE
unsetenv HELYX_BUILD_TYPE
unsetenv HELYX_PROJECT_VERSION
unsetenv HELYX_CC
unsetenv HELYX_CXX
unsetenv THIRDPARTY_INSTALL_DIR


#------------------------------------------------------------------------------
# Unset FOAM_* environment variables

unsetenv FOAM_CONFIG
unsetenv FOAM_CONFIG_ETC
unsetenv FOAM_CONFIG_MODE
unsetenv FOAM_SIGFPE
unsetenv FOAM_USER_APPBIN
unsetenv FOAM_USER_LIBBIN


#------------------------------------------------------------------------------
# Unset Third-Party environment variables

unsetenv HELYX_MPI_NAME
unsetenv MPI_ARCH_PATH
unsetenv MPI_BUFFER_SIZE
unsetenv OPAL_PREFIX
unsetenv I_MPI_ROOT
## Cleanup mpi prefix values if set to one of the paths on foamOldDirs
#if ( $?foamClean ) then
#    # openmpi:
#    if ( "`$foamClean -env=OPAL_PREFIX $foamOldDirs`" == "" ) unsetenv OPAL_PREFIX
#    # intelmpi:
#    if ( "`$foamClean -env=I_MPI_ROOT $foamOldDirs`" == "" ) unsetenv I_MPI_ROOT
#endif


#------------------------------------------------------------------------------
# Cleanup environment
# PATH, LD_LIBRARY_PATH, MANPATH

if ( $?foamClean ) then
    eval `$foamClean -csh-env=PATH "$foamOldDirs"`
    if ($?MANPATH) then
        eval `$foamClean -csh-env=MANPATH "$foamOldDirs"`
    endif
    if ($?LD_LIBRARY_PATH) then
        eval `$foamClean -csh-env=LD_LIBRARY_PATH "$foamOldDirs"`
    endif
endif

if ($?MANPATH) then
    if ("${MANPATH}" == "") unsetenv MANPATH
endif
if ($?LD_LIBRARY_PATH) then
    if ("${LD_LIBRARY_PATH}" == "") unsetenv LD_LIBRARY_PATH
endif


#------------------------------------------------------------------------------
# Cleanup aliases

unalias src
unalias run
unalias tut
unalias sol
unalias app
unalias util
unalias foam
unalias helyx


#------------------------------------------------------------------------------
# Cleanup auto-completions

# Remove old completions, which look like:
#   complete APPNAME 'p,*,`bash $HELYX_PROJECT_DIR/etc/ ...
if ($?prompt && $?tcsh) then  # Interactive tcsh only
    foreach cleaned (`complete | sed -ne '/HELYX_PROJECT/s/\t.*$//p'`)
        uncomplete $cleaned
    end
endif


#------------------------------------------------------------------------------
# Intermediate variables (do as last for a clean exit code)

unset cleaned foamClean foamOldDirs _foamFoundDir

#------------------------------------------------------------------------------
