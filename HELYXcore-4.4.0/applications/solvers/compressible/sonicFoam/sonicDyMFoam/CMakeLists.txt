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

set(TARGET_NAME "sonicDyMFoam-exe")

add_helyx_executable(${TARGET_NAME}
./sonicDyMFoam.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/applications/solvers/compressible/sonicFoam/sonicFoam
        ${HELYX_PROJECT_DIR}/applications/solvers/compressible/rhoPimpleFoam/rhoPimpleDyMFoam
        ${HELYX_PROJECT_DIR}/src/finiteVolume/cfdTools
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE sampling
    PRIVATE fvOptions
    PRIVATE helyxSolve
    PRIVATE dynamicFvMesh
    PRIVATE topoChangerFvMesh
    PRIVATE dynamicMesh
    PRIVATE meshTools
    PRIVATE OpenFOAM
)
