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

add_subdirectory(decompose)

add_subdirectory(decompositionMethods)

# N.B.  ${LIBRARY}DECOMP variables are being set on thirdPartyConfigurations.cmake.
# It looks like some libs linking against the *decomp libs are being compiled before
# src/parallel/decompose, so we need ${LIBRARY}DECOMP to be defined before we reach
# this file.
# Setting ${LIBRARY}DECOMP only here makes decomposition methods to not be
# recognized by some applications e.g., decomposePar.

if("${PTSCOTCH_FOUND}" AND "${MPI_FOUND}")
    add_subdirectory(ptscotchDecomp)
else()
    add_dummy_target(ptscotchDecomp "PTSCOTCH")
endif ()

if("${SCOTCH_FOUND}" AND "${MPI_FOUND}")
    add_subdirectory(scotchDecomp)
else()
    add_dummy_target(scotchDecomp "SCOTCH")
endif ()

if("${KAHIP_FOUND}" AND "${MPI_FOUND}")
    add_subdirectory(kahipDecomp)
else()
    add_dummy_target(kahipDecomp "KAHIP")
endif ()

if("${PARHIP_FOUND}" AND "${MPI_FOUND}")
    add_subdirectory(parhipDecomp)
else()
    add_dummy_target(parhipDecomp "PARHIP")
endif ()
