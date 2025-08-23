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
    (c) 2019-2020 Engys Ltd.

Description
    This file sets defaults for advanced and non-critical thirdParty variables.
    Note that no default paths are provided - all paths must be defined in
    userSettings.cmake (or equivalent).

[----------------------------------------------------------------------------]]

# Just the optional libraries
set(OPTIONAL_THIRDPARTY_LIBRARIES
# Dependencies of ThirdParties must be listed first
BOOST;  # Required by OPENVDB, PRECICE and sloanRenumber
CBLOSC; # Required by OPENVDB
TBB;    # Required by OPENVDB, (optionally) by OPENCASCADE, TBBbench, OpenFOAM, FiniteVolume
IMATH;  # Required by ALEMBIC
# ThirdParties
FFTW;
SCOTCH;
PTSCOTCH;
KAHIP;
PARHIP;
OPENVDB;
OPENCASCADE;
ALEMBIC;
PRECICE;
)

set(PARALLEL_THIRDPARTY_LIBRARIES
PTSCOTCH;
PARHIP;
PRECICE;
)

# Need to get the name of the MPI folder, which is used instead of the FOAM_MPI
# variable
get_filename_component(HELYX_MPI_NAME "${MPI_ARCH_PATH}" NAME CACHE)
