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

set(TARGET_NAME "eulerianMultiphaseSystems")

add_helyx_library(${TARGET_NAME} SHARED
./eulerianMultiphaseSystem/eulerianMultiphaseSystem.C
./eulerianMultiphaseSystem/newEulerianMultiphaseSystem.C
./eulerianMultiphaseSystem/eulerianMultiphaseSystems.C

./eulerianPhaseSystems/eulerianPhaseModel/eulerianPhaseModel/eulerianPhaseModel.C
./eulerianPhaseSystems/eulerianPhaseModel/eulerianPhaseModel/newEulerianPhaseModel.C
./eulerianPhaseSystems/eulerianPhaseModel/eulerianPhaseModel/eulerianPhaseModels.C
./eulerianPhaseSystems/eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.C
./eulerianPhaseSystems/eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.C
./eulerianPhaseSystems/eulerianPhasePair/orderedEulerianPhasePair/orderedEulerianPhasePair.C
./eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.C
./eulerianPhaseSystems/eulerianDiameterModels/eulerianDiameterModel/eulerianDiameterModel.C
./eulerianPhaseSystems/eulerianDiameterModels/eulerianDiameterModel/newEulerianDiameterModel.C
./eulerianPhaseSystems/eulerianDiameterModels/constantDiameter/constantDiameter.C
./eulerianPhaseSystems/eulerianDiameterModels/isothermalDiameter/isothermalDiameter.C
./eulerianPhaseSystems/EulerianBlendedInterfacialModel/eulerianBlendingMethods/eulerianBlendingMethod/eulerianBlendingMethod.C
./eulerianPhaseSystems/EulerianBlendedInterfacialModel/eulerianBlendingMethods/eulerianBlendingMethod/newEulerianBlendingMethod.C
./eulerianPhaseSystems/EulerianBlendedInterfacialModel/eulerianBlendingMethods/noBlending/noBlending.C
./eulerianPhaseSystems/EulerianBlendedInterfacialModel/eulerianBlendingMethods/linear/linear.C
./eulerianPhaseSystems/EulerianBlendedInterfacialModel/eulerianBlendingMethods/hyperbolic/hyperbolic.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/eulerianDragModel/eulerianDragModel.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/eulerianDragModel/newEulerianDragModel.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/Beetstra/Beetstra.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/segregated/segregated.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/Ergun/Ergun.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/Gibilaro/Gibilaro.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/Lain/Lain.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/SchillerNaumann/SchillerNaumann.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/SyamlalOBrien/SyamlalOBrien.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/TomiyamaCorrelated/TomiyamaCorrelated.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/TomiyamaAnalytic/TomiyamaAnalytic.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/Tenneti/Tenneti.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/TomiyamaKataokaZunSakaguchi/TomiyamaKataokaZunSakaguchi.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/WenYu/WenYu.C
./eulerianPhaseSystems/interfacialModels/eulerianDragModels/IshiiZuber/IshiiZuber.C
./eulerianPhaseSystems/interfacialModels/eulerianSwarmCorrections/eulerianSwarmCorrection/eulerianSwarmCorrection.C
./eulerianPhaseSystems/interfacialModels/eulerianSwarmCorrections/eulerianSwarmCorrection/newEulerianSwarmCorrection.C
./eulerianPhaseSystems/interfacialModels/eulerianSwarmCorrections/noSwarm/noSwarm.C
./eulerianPhaseSystems/interfacialModels/eulerianSwarmCorrections/TomiyamaSwarm/TomiyamaSwarm.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/eulerianLiftModel/eulerianLiftModel.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/eulerianLiftModel/newEulerianLiftModel.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/noLift/noLift.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/constantLiftCoefficient/constantLiftCoefficient.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/Moraga/Moraga.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/LegendreMagnaudet/LegendreMagnaudet.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/TomiyamaLift/TomiyamaLift.C
./eulerianPhaseSystems/interfacialModels/eulerianLiftModels/wallDampedLift/wallDampedLift.C
./eulerianPhaseSystems/interfacialModels/eulerianHeatTransferModels/eulerianHeatTransferModel/eulerianHeatTransferModel.C
./eulerianPhaseSystems/interfacialModels/eulerianHeatTransferModels/eulerianHeatTransferModel/eulerianHeatTransferModel.C
./eulerianPhaseSystems/interfacialModels/eulerianHeatTransferModels/eulerianHeatTransferModel/newEulerianHeatTransferModel.C
./eulerianPhaseSystems/interfacialModels/eulerianHeatTransferModels/RanzMarshall/RanzMarshall.C
./eulerianPhaseSystems/interfacialModels/eulerianHeatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.C
./eulerianPhaseSystems/interfacialModels/eulerianVirtualMassModels/eulerianVirtualMassModel/eulerianVirtualMassModel.C
./eulerianPhaseSystems/interfacialModels/eulerianVirtualMassModels/eulerianVirtualMassModel/newEulerianVirtualMassModel.C
./eulerianPhaseSystems/interfacialModels/eulerianVirtualMassModels/noVirtualMass/noVirtualMass.C
./eulerianPhaseSystems/interfacialModels/eulerianVirtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.C
./eulerianPhaseSystems/interfacialModels/eulerianVirtualMassModels/Lamb/Lamb.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/eulerianWallLubricationModel/eulerianWallLubricationModel.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/eulerianWallLubricationModel/newEulerianWallLubricationModel.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/noWallLubrication/noWallLubrication.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/Antal/Antal.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/Frank/Frank.C
./eulerianPhaseSystems/interfacialModels/eulerianWallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/eulerianTurbulentDispersionModel/eulerianTurbulentDispersionModel.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/eulerianTurbulentDispersionModel/newEulerianTurbulentDispersionModel.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/Burns/Burns.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/Gosman/Gosman.C
./eulerianPhaseSystems/interfacialModels/eulerianTurbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/eulerianAspectRatioModel/eulerianAspectRatioModel.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/eulerianAspectRatioModel/newEulerianAspectRatioModel.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/constantAspectRatio/constantAspectRatio.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.C
./eulerianPhaseSystems/interfacialModels/eulerianAspectRatioModels/Wellek/Wellek.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDependentModel/eulerianWallDependentModel.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/eulerianWallDampingModel/eulerianWallDampingModel.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/eulerianWallDampingModel/newEulerianWallDampingModel.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/noWallDamping/noWallDamping.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/interpolated/interpolatedWallDamping.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/linear/linearWallDamping.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/cosine/cosineWallDamping.C
./eulerianPhaseSystems/interfacialModels/eulerianWallDampingModels/sine/sineWallDamping.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianMassTransferModels/eulerianMassTransferModel/eulerianMassTransferModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianMassTransferModels/eulerianMassTransferModel/newEulerianMassTransferModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianMassTransferModels/Frossling/Frossling.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianMassTransferModels/sphericalMassTransfer/sphericalMassTransfer.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/eulerianInterfaceCompositionModel/eulerianInterfaceCompositionModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/eulerianInterfaceCompositionModel/newEulerianInterfaceCompositionModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/Henry/Henry.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/nonRandomTwoLiquid/nonRandomTwoLiquid.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/Raoult/Raoult.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianInterfaceCompositionModels/saturated/saturated.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/eulerianSaturationModel/eulerianSaturationModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/eulerianSaturationModel/newEulerianSaturationModel.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/Antoine/Antoine.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/AntoineExtended/AntoineExtended.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/ArdenBuck/ArdenBuck.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/polynomial/polynomial.C
./eulerianPhaseSystems/interfacialCompositionModels/eulerianSaturationModels/constantSaturationConditions/constantSaturationConditions.C

./multiphaseCompressibleTurbulenceModels/multiphaseCompressibleTurbulenceModels.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/transportModels
    PUBLIC
        turbulenceModels-includes
        turbulenceModelSchemes-includes
        combustionModels-includes
        chemistryModel-includes
        dynamicFvMesh-includes
        helyxSolve-includes
        phaseCompressibleTurbulenceModel-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fluidThermophysicalModels
    PRIVATE solidThermo
    PRIVATE specie
    PUBLIC thermophysicalMixtures
    PUBLIC materialModels
    PUBLIC multiphaseMaterials
    PRIVATE thermophysicalProperties
    PRIVATE compressibleTurbulenceModels
    PRIVATE combustionModels
    PRIVATE finiteVolume
    PRIVATE finiteArea
    PRIVATE meshTools
)
