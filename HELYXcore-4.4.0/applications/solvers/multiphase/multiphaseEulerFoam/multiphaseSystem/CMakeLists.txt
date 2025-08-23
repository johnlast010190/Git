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

set(TARGET_NAME "multiphaseSystem")

add_helyx_library(${TARGET_NAME} SHARED
./phaseModel/phaseModel.C
./diameterModels/diameterModel/diameterModel.C
./diameterModels/diameterModel/newDiameterModel.C
./diameterModels/constantDiameter/constantDiameter.C
./diameterModels/isothermalDiameter/isothermalDiameter.C
./interfacialModels/dragModels/dragModel/dragModel.C
./interfacialModels/dragModels/dragModel/newDragModel.C
./interfacialModels/dragModels/Ergun/Ergun.C
./interfacialModels/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
./interfacialModels/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
./interfacialModels/dragModels/SchillerNaumann/SchillerNaumann.C
./interfacialModels/dragModels/Gibilaro/Gibilaro.C
./interfacialModels/dragModels/WenYu/WenYu.C
./interfacialModels/dragModels/SyamlalOBrien/SyamlalOBrien.C
./interfacialModels/dragModels/blended/blended.C
./interfacialModels/dragModels/interface/interface.C
./interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.C
./interfacialModels/heatTransferModels/heatTransferModel/newHeatTransferModel.C
./interfacialModels/heatTransferModels/RanzMarshall/RanzMarshall.C
./alphaContactAngle/alphaContactAngleMPFvPatchScalarField.C
./multiphaseSystem.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/applications/solvers/multiphase/multiphaseEulerFoam/multiphaseSystem/alphaContactAngle
        ${HELYX_PROJECT_DIR}/src/transportModels
        incompressibleTransportModels-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE interfaceProperties
    PRIVATE incompressibleTransportModels
    PRIVATE finiteVolume
    PRIVATE meshTools
)
