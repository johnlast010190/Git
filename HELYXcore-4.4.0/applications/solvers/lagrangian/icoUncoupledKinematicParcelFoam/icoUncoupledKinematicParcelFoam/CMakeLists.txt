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
    (c) 2019-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "icoUncoupledKinematicParcelFoam-exe")

add_helyx_executable(${TARGET_NAME}
./icoUncoupledKinematicParcelFoam.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fluidThermophysicalModels-includes
        ${HELYX_PROJECT_DIR}/src/transportModels
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE lagrangian
    PRIVATE lagrangianIntermediate
    PRIVATE lagrangianTurbulence
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE radiationModels
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE incompressibleTransportModels
    PRIVATE finiteVolume
    PRIVATE fvOptions
    PRIVATE helyxSolve
    PRIVATE meshTools
    PRIVATE regionModels
    PRIVATE OpenFOAM
)
