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

set(TARGET_NAME "fvOptions")

add_helyx_library(${TARGET_NAME} SHARED
./cellSetOption/cellSetOption.C
./cellSetOption/cellSetOptionIO.C
./cellSetOption/sweptCellMRF/sweptCellMRF.C
./interRegionOption/interRegionOption.C
./interRegionOption/interRegionOptionIO.C
./sources/general/codedSource/codedSource.C
./sources/general/semiImplicitSource/semiImplicitSource.C
./sources/general/timeDependentBoundarySource/timeDependentBoundarySource.C
./sources/derived/accelerationSource/accelerationSource.C
./sources/derived/acousticDampingSource/acousticDampingSource.C
./sources/derived/acousticSource/acousticSource.C
./sources/derived/actuationDiskSource/actuationDiskSource.C
./sources/derived/buoyancyForce/buoyancyForce.C
./sources/derived/buoyancyForce/buoyancyForceIO.C
./sources/derived/buoyancyEnergy/buoyancyEnergy.C
./sources/derived/buoyancyEnergy/buoyancyEnergyIO.C
./sources/derived/directionalPressureGradientExplicitSource/directionalPressureGradientExplicitSource.C
./sources/derived/directionalPressureGradientExplicitSource/directionalPressureGradientExplicitSourceIO.C
./sources/derived/effectivenessHeatExchangerSource/effectivenessHeatExchangerSource.C
./sources/derived/externalHeatExchangerSource/externalHeatExchangerSource.C
./sources/derived/explicitPorositySource/explicitPorositySource.C
./sources/derived/frameSources/frameSources.C
./sources/derived/MRFSource/MRFSource.C
./sources/derived/MRFSource/MRFSourceList.C
./sources/derived/MRFSource/IOMRFSourceList.C
./sources/derived/meanVelocityForce/meanVelocityForce.C
./sources/derived/meanVelocityForce/meanVelocityForceIO.C
./sources/derived/meanVelocityForce/patchMeanVelocityForce/patchMeanVelocityForce.C
./sources/derived/radialActuationDiskSource/radialActuationDiskSource.C
./sources/derived/rotatingVelocityForce/rotatingVelocityForceIO.C
./sources/derived/rotatingVelocityForce/rotatingVelocityForce.C
./sources/derived/rotatingVelocityForce/patchRotatingVelocityForce/patchRotatingVelocityForce.C
./sources/derived/rotorDiskSource/rotorDiskSource.C
./sources/derived/rotorDiskSource/bladeModel/bladeModel.C
./sources/derived/rotorDiskSource/profileModel/profileModel.C
./sources/derived/rotorDiskSource/profileModel/profileModelList.C
./sources/derived/rotorDiskSource/profileModel/lookup/lookupProfile.C
./sources/derived/rotorDiskSource/profileModel/series/seriesProfile.C
./sources/derived/rotorDiskSource/trimModel/trimModel/trimModel.C
./sources/derived/rotorDiskSource/trimModel/trimModel/trimModelNew.C
./sources/derived/rotorDiskSource/trimModel/fixed/fixedTrim.C
./sources/derived/rotorDiskSource/trimModel/targetCoeff/targetCoeffTrim.C
./sources/derived/solidificationMeltingSource/solidificationMeltingSource.C
./sources/derived/solidificationMeltingSource/solidificationMeltingSourceIO.C
./sources/derived/tabulatedAccelerationSource/tabulatedAccelerationSource.C
./sources/derived/tabulatedAccelerationSource/tabulated6DoFAcceleration/tabulated6DoFAcceleration.C
./sources/derived/thermalSource/thermalSource.C
./sources/derived/PatersonAnnularThrusterSource/PatersonAnnularThrusterSource.C
./sources/derived/thrusterActuationDiskSource/thrusterActuationDiskSource.C
./sources/derived/constantThrustActuationDiskSource/constantThrustActuationDiskSource.C
./sources/derived/radialThrustActuationDiskSource/controllers/controllers.C
./sources/derived/radialThrustActuationDiskSource/radialThrustActuationDiskSource.C
./sources/derived/forceActuationDiskSource/forceActuationDiskSource.C
./sources/derived/exponentialThermalSource/exponentialThermalSource.C
./sources/derived/velocityLimiterSource/velocityLimiterSource.C
./sources/derived/scalarFieldLimiterSource/scalarFieldLimiterSource.C
./sources/derived/artificialViscosityUSource/artificialViscosityUSource.C
./sources/derived/GRFSource/GRFSource.C
./sources/derived/GRFSource/GRF/GRF.C
./sources/derived/solidRotationSource/solidRotationSource.C
./sources/derived/vegetationSource/vegetationSource.C
./sources/derived/vegetationSource/dalpeMassonCanopySource/dalpeMassonCanopySource.C
./sources/derived/vegetationSource/svenssonHaggkvistCanopySource/svenssonHaggkvistCanopySource.C
./sources/derived/absorptionEmissionSource/absorptionEmissionSource.C
./sources/derived/porousAbsorptionEmissionSource/porousAbsorptionEmissionSource.C
./sources/derived/thermalRunawaySource/thermalRunawaySource.C
./sources/derived/viscousDissipation/viscousDissipation.C
./sources/derived/settlingVelocityTakacs/settlingVelocityTakacs.C
./sources/derived/turbulenceControlledDecaySource/turbulenceControlledDecaySource.C
./sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModel.C
./sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModelIO.C
./sources/interRegion/interRegionHeatTransfer/constantHeatTransfer/constantHeatTransfer.C
./sources/interRegion/interRegionHeatTransfer/tabulatedHeatTransfer/tabulatedHeatTransfer.C
./sources/interRegion/interRegionHeatTransfer/tabulatedNTUHeatTransfer/tabulatedNTUHeatTransfer.C
./sources/interRegion/interRegionHeatTransfer/variableHeatTransfer/variableHeatTransfer.C
./sources/interRegion/interRegionHeatTransfer/tabulatedEffectivenessHeatTransfer/tabulatedEffectivenessHeatTransfer.C
./sources/interRegion/interRegionHeatTransfer/tabulatedUALHeatTransfer/tabulatedUALHeatTransfer.C
./sources/interRegion/interRegionExplicitPorositySource/interRegionExplicitPorositySource.C
./constraints/general/fixedValueConstraint/fixedValueConstraints.C
./constraints/general/fieldNeighbourValueConstraint/fieldNeighbourValueConstraints.C
./constraints/general/explicitSetDataEntryValue/explicitSetDataEntryValue.C
./constraints/derived/fixedTemperatureConstraint/fixedTemperatureConstraint.C
./constraints/derived/velocityDampingConstraint/velocityDampingConstraint.C
./constraints/derived/temperatureLimitsConstraint/temperatureLimitsConstraint.C
./transformations/anisotropicConductivity/anisotropicConductivityTransformation.C
./corrections/limitTemperature/limitTemperature.C
./corrections/limitVelocity/limitVelocity.C
./corrections/CourantNumberLimit/CourantNumberLimit.C
./corrections/atcMask/atcMask.C
./postProcessing/equationFlux/equationFluxes.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE sampling
    PRIVATE meshTools
    PRIVATE dynamicFvMesh
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE forces
    PRIVATE radiationModels
    PRIVATE incompressibleTransportModels
    PRIVATE eulerianMultiphaseSystems
    PRIVATE multiphaseMaterials
    PRIVATE materialModels
    PRIVATE solidThermo
)
