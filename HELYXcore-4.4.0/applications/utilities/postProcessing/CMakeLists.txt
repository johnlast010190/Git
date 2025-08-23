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
    (c) 2020-2024 Engys Ltd.

[----------------------------------------------------------------------------]]

add_subdirectory(dataConversion)

add_subdirectory(graphics)

add_subdirectory(lagrangian)

add_subdirectory(helyxSample)

add_subdirectory(levelSet)

add_subdirectory(mkSymmetric)

add_subdirectory(np2foam)

add_subdirectory(objectMask)

add_subdirectory(sampleGeom)

add_subdirectory(mapObject)

add_subdirectory(mapGeom)

add_subdirectory(miscellaneous)

add_subdirectory(postProcess)

add_subdirectory(scalarField)

add_subdirectory(wall)

if(${FFTW_FOUND})
    add_subdirectory(noise)
    add_subdirectory(volumetricFFT)
else()
    add_dummy_target(noise-exe FFTW)
    add_dummy_target(volumetricFFT-exe FFTW)
endif()



