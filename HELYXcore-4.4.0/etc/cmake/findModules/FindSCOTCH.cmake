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
    Custom findModule for scotch

[----------------------------------------------------------------------------]]


set(SCOTCH_PATH_SUFFIXES "../lib;${SCOTCH_PATH_SUFFIXES}")
helyx_find_thirdparty_package(SCOTCH "scotch;scotcherrexit" scotch.h)

# If OpenCASCADE compilation is disabled, simply return from this script.
# There is no need to evaluate dependencies
if(SCOTCH_REQUIRED STREQUAL "OFF" OR NOT SCOTCH_FOUND)
    return()
endif()

# Scotch has an implicit dependency on zlib, so connect the imported targets together.
# thirdparty_blah is just a list of imported target names. Yes, this works. No, it doesn't call the linker :D.
target_link_libraries(scotch INTERFACE "${THIRDPARTY_ZLIB}")
