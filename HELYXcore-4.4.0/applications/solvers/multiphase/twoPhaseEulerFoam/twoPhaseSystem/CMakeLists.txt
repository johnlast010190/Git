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

set(TARGET_NAME "compressibleTwoPhaseSystem")

add_helyx_library(${TARGET_NAME} SHARED
./phaseModel/phaseModel.C
./diameterModels/diameterModel/diameterModel.C
./diameterModels/diameterModel/newDiameterModel.C
./diameterModels/constantDiameter/constantDiameter.C
./diameterModels/isothermalDiameter/isothermalDiameter.C
./diameterModels/IATE/IATE.C
./diameterModels/IATE/IATEsources/IATEsource/IATEsource.C
./diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.C
./diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.C
./diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.C
./BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethod.C
./BlendedInterfacialModel/blendingMethods/blendingMethod/newBlendingMethod.C
./BlendedInterfacialModel/blendingMethods/noBlending/noBlending.C
./BlendedInterfacialModel/blendingMethods/linear/linear.C
./BlendedInterfacialModel/blendingMethods/hyperbolic/hyperbolic.C
./phasePair/phasePair/phasePair.C
./phasePair/orderedPhasePair/orderedPhasePair.C
./twoPhaseSystem.C
./phaseCompressibleTurbulenceModel.C
./interfacialModels/dragModels/dragModel/dragModel.C
./interfacialModels/dragModels/dragModel/newDragModel.C
./interfacialModels/dragModels/segregated/segregated.C
./interfacialModels/dragModels/Ergun/Ergun.C
./interfacialModels/dragModels/Gibilaro/Gibilaro.C
./interfacialModels/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
./interfacialModels/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
./interfacialModels/dragModels/Lain/Lain.C
./interfacialModels/dragModels/SchillerNaumann/SchillerNaumann.C
./interfacialModels/dragModels/SyamlalOBrien/SyamlalOBrien.C
./interfacialModels/dragModels/TomiyamaCorrelated/TomiyamaCorrelated.C
./interfacialModels/dragModels/TomiyamaAnalytic/TomiyamaAnalytic.C
./interfacialModels/dragModels/WenYu/WenYu.C
./interfacialModels/dragModels/IshiiZuber/IshiiZuber.C
./interfacialModels/swarmCorrections/swarmCorrection/swarmCorrection.C
./interfacialModels/swarmCorrections/swarmCorrection/newSwarmCorrection.C
./interfacialModels/swarmCorrections/noSwarm/noSwarm.C
./interfacialModels/swarmCorrections/TomiyamaSwarm/TomiyamaSwarm.C
./interfacialModels/liftModels/liftModel/liftModel.C
./interfacialModels/liftModels/liftModel/newLiftModel.C
./interfacialModels/liftModels/noLift/noLift.C
./interfacialModels/liftModels/constantLiftCoefficient/constantLiftCoefficient.C
./interfacialModels/liftModels/Moraga/Moraga.C
./interfacialModels/liftModels/LegendreMagnaudet/LegendreMagnaudet.C
./interfacialModels/liftModels/TomiyamaLift/TomiyamaLift.C
./interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.C
./interfacialModels/heatTransferModels/heatTransferModel/newHeatTransferModel.C
./interfacialModels/heatTransferModels/RanzMarshall/RanzMarshall.C
./interfacialModels/heatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.C
./interfacialModels/virtualMassModels/virtualMassModel/virtualMassModel.C
./interfacialModels/virtualMassModels/virtualMassModel/newVirtualMassModel.C
./interfacialModels/virtualMassModels/noVirtualMass/noVirtualMass.C
./interfacialModels/virtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.C
./interfacialModels/virtualMassModels/Lamb/Lamb.C
./interfacialModels/wallLubricationModels/wallLubricationModel/wallLubricationModel.C
./interfacialModels/wallLubricationModels/wallLubricationModel/newWallLubricationModel.C
./interfacialModels/wallLubricationModels/noWallLubrication/noWallLubrication.C
./interfacialModels/wallLubricationModels/Antal/Antal.C
./interfacialModels/wallLubricationModels/Frank/Frank.C
./interfacialModels/wallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.C
./interfacialModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.C
./interfacialModels/turbulentDispersionModels/turbulentDispersionModel/newTurbulentDispersionModel.C
./interfacialModels/turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.C
./interfacialModels/turbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.C
./interfacialModels/turbulentDispersionModels/Burns/Burns.C
./interfacialModels/turbulentDispersionModels/Gosman/Gosman.C
./interfacialModels/turbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.C
./interfacialModels/aspectRatioModels/aspectRatioModel/aspectRatioModel.C
./interfacialModels/aspectRatioModels/aspectRatioModel/newAspectRatioModel.C
./interfacialModels/aspectRatioModels/constantAspectRatio/constantAspectRatio.C
./interfacialModels/aspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.C
./interfacialModels/aspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.C
./interfacialModels/aspectRatioModels/Wellek/Wellek.C
./interfacialModels/wallDependentModel/wallDependentModel.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/applications/solvers/multiphase/twoPhaseEulerFoam/twoPhaseSystem
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/transportModel
        phaseCompressibleTurbulenceModel-includes
        meshTools-includes
        sampling-includes
        phaseCompressibleTurbulenceModel-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE incompressibleTransportModels
    PRIVATE fluidThermophysicalModels
    PRIVATE specie
    PRIVATE finiteVolume
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE fvOptions
    PRIVATE OpenFOAM
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/BlendedInterfacialModel/blendingMethods/linear/linear.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
