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
    Custom findModule for ptscotch

[----------------------------------------------------------------------------]]


# N.B.  If using system scotch and system MPI, ptscotch may be installed in the
# system MPI directory.  helyx_find_thirdparty_package() handles this, but it
# should perhaps be handled here instead.
set(PTSCOTCH_PATH_SUFFIXES "../lib/${HELYX_MPI_NAME};../lib;../${HELYX_MPI_NAME};${PTSCOTCH_PATH_SUFFIXES}")
helyx_find_thirdparty_package(PTSCOTCH "ptscotch;ptscotcherrexit" ptscotch.h)
