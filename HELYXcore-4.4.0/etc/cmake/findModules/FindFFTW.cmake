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
    Custom findModule for FFTW

[----------------------------------------------------------------------------]]


# No need to do this, because FFTW doesn't have a find module!
# Push CMAKE_MODULE_PATH to enable calling CMake's find_package() from here
# set(STORED_CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}")
# set(CMAKE_MODULE_PATH "")

# We should reject stactic links of fftw, so
#forcing search for the dynamic library
if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    set(fftw_dynamic_lib "libfftw3.dll")
else()
    set(fftw_dynamic_lib "libfftw3.so")
endif()

helyx_find_thirdparty_package(FFTW ${fftw_dynamic_lib} fftw3.h)

# No need to do this, because FFTW doesn't have a find module!
# Pop CMAKE_MODULE_PATH
# set(CMAKE_MODULE_PATH "${STORED_CMAKE_MODULE_PATH}")
