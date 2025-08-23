#!/bin/bash
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
#    (c) 2021 Engys Ltd.
#
#Description
#     Autocompletion information for emake.
#------------------------------------------------------------------------------


FLAGS="
-refreshCache
-noBuild
-cleanFirst
-clean
-thirdParty
-availableTargets
-generateSourceableFiles
-help
--refactor-cpp-includes
--refactor-cmake-includes
--refresh-cache
--no-build
--clean-first
--clean
--third-party
--available-targets
--generate-sourceable-files
--submodules
--pull
--help
"

THIRDPARTY_TARGETS="
thirdParty_MPI
thirdParty_SCOTCH
thirdParty_KAHIP
thirdParty_PARHIP
thirdParty_FFTW
thirdParty_TBB
thirdParty_CBLOSC
thirdParty_BOOST
thirdParty_OPENVDB
thirdParty_OPENCASCADE
thirdParty_IMATH
thirdParty_ALEMBIC
thirdParty_PRECICE
"

TARGETS="
all
install
package-binaries
packageBinaries
package-source
packageSource
@TARGETS_STRING@
"

_emake()
{
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    if [[ ${cur} == -* ]] ; then
        COMPREPLY=( $(compgen -W "${FLAGS}" -- ${cur}) )
    elif [[ $prev == -t ]]; then
        COMPREPLY=( $(compgen -W "$THIRDPARTY_TARGETS" -- ${cur}) )
    else
        COMPREPLY=( $(compgen -W "${TARGETS}" -- ${cur}) )
    fi
    return 0
}

complete -F _emake emake
