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

set(TARGET_NAME "lagrangianIntermediate")

add_helyx_library(${TARGET_NAME} SHARED
./submodels/CloudFunctionObjects/CloudToVTK/vtkTools.C
./clouds/baseClasses/kinematicCloud/kinematicCloud.C
./clouds/baseClasses/thermoCloud/thermoCloud.C
./clouds/baseClasses/reactingCloud/reactingCloud.C
./clouds/baseClasses/reactingMultiphaseCloud/reactingMultiphaseCloud.C
./parcels/derived/basicKinematicParcel/defineBasicKinematicParcel.C
./parcels/derived/basicKinematicParcel/makeBasicKinematicParcelSubmodels.C
./parcels/derived/basicKinematicCollidingParcel/defineBasicKinematicCollidingParcel.C
./parcels/derived/basicKinematicCollidingParcel/makeBasicKinematicCollidingParcelSubmodels.C
./parcels/derived/basicThermoParcel/defineBasicThermoParcel.C
./parcels/derived/basicThermoParcel/makeBasicThermoParcelSubmodels.C
./parcels/derived/basicReactingParcel/defineBasicReactingParcel.C
./parcels/derived/basicReactingParcel/makeBasicReactingParcelSubmodels.C
./parcels/derived/basicReactingMultiphaseParcel/defineBasicReactingMultiphaseParcel.C
./parcels/derived/basicReactingMultiphaseParcel/makeBasicReactingMultiphaseParcelSubmodels.C
./parcels/derived/basicKinematicMPPICParcel/defineBasicKinematicMPPICParcel.C
./parcels/derived/basicKinematicMPPICParcel/makeBasicKinematicMPPICParcelSubmodels.C
./submodels/addOns/radiation/absorptionEmission/cloudAbsorptionEmission/cloudAbsorptionEmission.C
./submodels/addOns/radiation/scatter/cloudScatter/cloudScatter.C
./submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.C
./submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionDataList.C
./submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionData.C
./submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionDataIO.C
./submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionDataIOList.C
./submodels/Kinematic/InjectionModel/PatchInjection/patchInjectionBase.C
./submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionData.C
./submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionDataIO.C
./submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionDataIOList.C
./submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionData.C
./submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionDataIO.C
./submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionDataIOList.C
./submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionData.C
./submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIO.C
./submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIOList.C
./submodels/MPPIC/ParticleStressModels/ParticleStressModel/ParticleStressModel.C
./submodels/MPPIC/ParticleStressModels/HarrisCrighton/HarrisCrighton.C
./submodels/MPPIC/ParticleStressModels/Lun/Lun.C
./submodels/MPPIC/ParticleStressModels/exponential/exponential.C
./submodels/MPPIC/CorrectionLimitingMethods/CorrectionLimitingMethod/CorrectionLimitingMethod.C
./submodels/MPPIC/CorrectionLimitingMethods/noCorrectionLimiting/noCorrectionLimiting.C
./submodels/MPPIC/CorrectionLimitingMethods/absolute/absolute.C
./submodels/MPPIC/CorrectionLimitingMethods/relative/relative.C
./submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.C
./submodels/MPPIC/TimeScaleModels/equilibrium/equilibrium.C
./submodels/MPPIC/TimeScaleModels/nonEquilibrium/nonEquilibrium.C
./submodels/MPPIC/TimeScaleModels/isotropic/isotropic.C
./IntegrationScheme/makeIntegrationSchemes.C
./phaseProperties/phaseProperties/phaseProperties.C
./phaseProperties/phaseProperties/phasePropertiesIO.C
./phaseProperties/phasePropertiesList/phasePropertiesList.C
./clouds/Templates/KinematicCloud/cloudSolution/cloudSolution.C
./submodels/MPPIC/AveragingMethods/makeAveragingMethods.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/transportModels
        surfMesh-includes
    PUBLIC
        regionModels-includes
        materialModels-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE lagrangian
    PRIVATE materialModels
    PRIVATE distributionModels
    PRIVATE specie
    PRIVATE fluidThermophysicalModels
    PRIVATE thermophysicalProperties
    PRIVATE SLGThermo
    PRIVATE radiationModels
    PRIVATE incompressibleTransportModels
    PRIVATE regionModels
    PRIVATE dynamicMesh
    PRIVATE dynamicFvMesh
    PRIVATE sampling
    PRIVATE finiteVolume
    PRIVATE meshTools
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicKinematicCollidingParcel/makeBasicKinematicCollidingParcelSubmodels.C
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicKinematicMPPICParcel/makeBasicKinematicMPPICParcelSubmodels.C
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicKinematicParcel/makeBasicKinematicParcelSubmodels.C
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicReactingMultiphaseParcel/makeBasicReactingMultiphaseParcelSubmodels.C
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicReactingParcel/makeBasicReactingParcelSubmodels.C
    ${CMAKE_CURRENT_LIST_DIR}/parcels/derived/basicThermoParcel/makeBasicThermoParcelSubmodels.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
