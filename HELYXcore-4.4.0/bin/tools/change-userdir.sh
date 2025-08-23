#---------------------------------------------------------------------------
#|       o        |
#|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
#|   o   O   o    |  Version : 4.4.0
#|    o     o     |  ENGYS Ltd. <http://engys.com/>
#|       o        |
#---------------------------------------------------------------------------
#License
#    This file is part of HELYXcore.
#    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
#
#    HELYXcore is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

#Copyright
#    (c) 2017 OpenCFD Ltd. 
#
# Script
#     . change-userdir.sh PREFIX [SUFFIX]
#
#     Shortcuts (prefix)
#         -home           "$HOME/OpenFOAM/$USER-$HELYX_PROJECT_VERSION"
#         -none           remove from environment
#
#     Shortcuts (suffix)
#         -platforms      "platforms/$HELYX_OPTIONS"
#
# Description
#     Change CMAKE_PROJECT_NAME_USER_DIR, HELYX_USER_APPBIN, HELYX_USER_LIBBIN
#     and the respective entries in PATH, LD_LIBRARY_PATH.
#     Also adjusts FOAM_RUN.
#
#     This can be useful with compiling additional OpenFOAM programs
#     (that use HELYX_USER_APPBIN, HELYX_USER_LIBBIN for their build),
#     to avoid conflicts with the normal user bin/lib files.
#
#     The suffix value should normally include "platforms/$HELYX_OPTIONS"
#
# Example
#     . /path/change-userdir.sh -home -platforms
#
#   corresponds to the standard user location:
#
#     $HOME/OpenFOAM/$USER-$HELYX_PROJECT_VERSION/platforms/$HELYX_OPTIONS
#
#------------------------------------------------------------------------------

if [ "$#" -ge 1 ]
then
    prefix="$1"
    suffix="$2"

    foamOldDirs="$HELYX_USER_APPBIN $HELYX_USER_LIBBIN"
    foamClean=$HELYX_PROJECT_DIR/bin/foamCleanPath
    if [ -x "$foamClean" ]
    then
        cleaned=$($foamClean "$PATH" "$foamOldDirs") && PATH="$cleaned"
        cleaned=$($foamClean "$LD_LIBRARY_PATH" "$foamOldDirs") \
            && LD_LIBRARY_PATH="$cleaned"
    fi

    case "$suffix" in
        -plat*) suffix="platforms/$HELYX_OPTIONS" ;;
    esac
    case "$prefix" in
        -home) prefix="$HOME/OpenFOAM/$USER-${HELYX_PROJECT_VERSION:-unknown}" ;;
        -none) unset prefix ;;
    esac

    if [ -n "$prefix" ]
    then
        export CMAKE_PROJECT_NAME_USER_DIR="$prefix"
        export FOAM_RUN="$prefix/run"

        prefix="$prefix${suffix:+/}${suffix}"
        export HELYX_USER_APPBIN="$prefix/bin"
        export HELYX_USER_LIBBIN="$prefix/lib"

        PATH="$HELYX_USER_APPBIN:$PATH"
        LD_LIBRARY_PATH="$HELYX_USER_LIBBIN:$LD_LIBRARY_PATH"
    else
        unset CMAKE_PROJECT_NAME_USER_DIR FOAM_RUN HELYX_USER_APPBIN HELYX_USER_LIBBIN
    fi
fi

unset foamClean foamOldDirs cleaned prefix suffix

#------------------------------------------------------------------------------
