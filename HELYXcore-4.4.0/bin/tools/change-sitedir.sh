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
#     . change-sitedir.sh PREFIX [SUFFIX]
#
#     Shortcuts (prefix)
#         -prefix         "$CMAKE_PROJECT_NAME_INST_DIR/site"
#         -project        "$HELYX_PROJECT_DIR/site"
#         -none           remove from environment
#
#     Shortcuts (suffix)
#         -platforms      "platforms/$HELYX_OPTIONS"
#
# Description
#     Change HELYX_CUSTOM_DATA_DIR, FOAM_SITE_APPBIN, FOAM_SITE_LIBBIN
#     and the respective entries in PATH, LD_LIBRARY_PATH.
#
#     This can be useful when temporarily reassigning the site directory
#     when packaging OpenFOAM.
#
#     The suffix value should normally include "platforms/$HELYX_OPTIONS"
#
# Example
#     . /path/change-sitedir.sh -prefix -platforms
#
#   corresponds to the standard site location:
#
#     $CMAKE_PROJECT_NAME_INST_DIR/site{/$HELYX_PROJECT_VERSION/platforms/$HELYX_OPTIONS}
#
#------------------------------------------------------------------------------

if [ "$#" -ge 1 ]
then
    prefix="$1"
    suffix="$2"

    foamOldDirs="$FOAM_SITE_APPBIN $FOAM_SITE_LIBBIN \
        $HELYX_CUSTOM_DATA_DIR $CMAKE_PROJECT_NAME_INST_DIR/site $HELYX_PROJECT_DIR/site"
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
        -prefix)  prefix="$CMAKE_PROJECT_NAME_INST_DIR/site" ;;
        -project) prefix="$HELYX_PROJECT_DIR/site" ;;
        -none)    unset prefix ;;
    esac

    if [ -n "$prefix" ]
    then
        export HELYX_CUSTOM_DATA_DIR="$prefix"

        prefix="$prefix/${HELYX_PROJECT_VERSION:-unknown}${suffix:+/}${suffix}"

        export FOAM_SITE_APPBIN="$prefix/bin"
        export FOAM_SITE_LIBBIN="$prefix/lib"
        PATH="$FOAM_SITE_APPBIN:$PATH"
        LD_LIBRARY_PATH="$FOAM_SITE_LIBBIN:$LD_LIBRARY_PATH"
    else
        unset HELYX_CUSTOM_DATA_DIR FOAM_SITE_APPBIN FOAM_SITE_LIBBIN
    fi
fi

unset foamClean foamOldDirs cleaned prefix suffix

#------------------------------------------------------------------------------
