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
    (c) 2019 Engys Ltd.

[----------------------------------------------------------------------------]]

add_subdirectory(basic)

add_subdirectory(combustion)

add_subdirectory(compressible)

add_subdirectory(discreteMethods)

if(${FFTW_FOUND})
    add_subdirectory(DNS)
else()
    add_dummy_target(dnsFoam-exe FFTW)
endif()

add_subdirectory(electromagnetics)

add_subdirectory(heatTransfer)

add_subdirectory(helyxSolve)

add_subdirectory(incompressible)

add_subdirectory(lagrangian)

add_subdirectory(multiphase)

add_subdirectory(stressAnalysis)

add_subdirectory(vibroacoustics)
