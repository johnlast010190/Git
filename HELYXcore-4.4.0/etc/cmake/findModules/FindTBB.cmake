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
    (c) 2022 Engys Ltd.

Description
    Custom findModule for TBB

[----------------------------------------------------------------------------]]


# N.B.  TBB headers are not in standard include dir, instead, they are located
# in */include/tbb/, but TBB own code expects includes from */include/.
# Inside 'helyx_find_thirdparty_package' the found includes path is corrected to
# point to */include/ (instead of */include/tbb/) and then, the headers are
# associated to the TBB libs.
# This is done by listing the headers in 'helyx_find_thirdparty_package' as a
# suffix of the expect include folder.
helyx_find_thirdparty_package(TBB "tbb;tbbmalloc" "tbb/tbb.h")
