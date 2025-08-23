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

add_subdirectory(MPPICInterFoam)

add_subdirectory(algebraicSlipFoam)

add_subdirectory(cavitatingFoam)

add_subdirectory(compressibleInterFoam)

add_subdirectory(compressibleMultiphaseInterFoam)

add_subdirectory(driftFluxFoam)

add_subdirectory(interFoam)

add_subdirectory(interIsoFoam)

add_subdirectory(interPhaseChangeFoam)

add_subdirectory(meltFoam)

add_subdirectory(multiphaseEulerFoam)

add_subdirectory(multiphaseInterFoam)

add_subdirectory(potentialFreeSurfaceFoam)

add_subdirectory(reactingEulerFoam)

add_subdirectory(twoLiquidMixingFoam)

add_subdirectory(twoPhaseEulerFoam)

add_subdirectory(multiphaseDriftFluxFoam)
