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

set(TARGET_NAME "PDRFoam-exe")

add_helyx_executable(${TARGET_NAME}
./XiModels/XiModel/XiModel.C
./XiModels/XiModel/XiModelNew.C
./XiModels/fixed/fixed.C
./XiModels/algebraic/algebraic.C
./XiModels/transport/transport.C
./XiModels/XiEqModels/XiEqModel/XiEqModel.C
./XiModels/XiEqModels/XiEqModel/XiEqModelNew.C
./XiModels/XiEqModels/Gulder/Gulder.C
./XiModels/XiEqModels/instabilityXiEq/instabilityXiEq.C
./XiModels/XiEqModels/SCOPEBlendXiEq/SCOPEBlendXiEq.C
./XiModels/XiEqModels/SCOPEXiEq/SCOPEXiEq.C
./XiModels/XiGModels/XiGModel/XiGModel.C
./XiModels/XiGModels/XiGModel/XiGModelNew.C
./XiModels/XiGModels/KTS/KTS.C
./XiModels/XiGModels/instabilityG/instabilityG.C
./PDRModels/turbulence/PDRkEpsilon/PDRkEpsilon.C
./PDRModels/dragModels/PDRDragModel/PDRDragModel.C
./PDRModels/dragModels/PDRDragModel/PDRDragModelNew.C
./PDRModels/dragModels/basic/basic.C
./PDRModels/XiEqModels/basicXiSubXiEq/basicXiSubXiEq.C
./PDRModels/XiGModels/basicXiSubG/basicXiSubG.C
./laminarFlameSpeed/SCOPE/SCOPELaminarFlameSpeed.C
./PDRFoam.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/applications/solvers/combustion/PDRFoam/XiModels
        ${HELYX_PROJECT_DIR}/applications/solvers/combustion/PDRFoam/XiModels/XiEqModels
        ${HELYX_PROJECT_DIR}/applications/solvers/combustion/PDRFoam/XiModels/XiGModels
        ${HELYX_PROJECT_DIR}/applications/solvers/combustion/PDRFoam/PDRModels/dragModels
        ${HELYX_PROJECT_DIR}/applications/solvers/combustion/PDRFoam/laminarFlameSpeed
        dynamicMesh-includes
        surfMesh-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE engine
    PRIVATE meshTools
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE laminarFlameSpeedModels
    PRIVATE finiteVolume
    PRIVATE dynamicFvMesh
    PRIVATE fvOptions
    PRIVATE helyxSolve
    PRIVATE OpenFOAM
)


set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/PDRModels/XiGModels/basicXiSubG/basicXiSubG.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
