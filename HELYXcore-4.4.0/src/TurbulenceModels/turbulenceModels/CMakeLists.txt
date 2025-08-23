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

set(TARGET_NAME "turbulenceModels")

add_helyx_library(${TARGET_NAME} SHARED
./turbulenceModel.C
./RAS/v2f/v2fBase.C
./LES/LESdeltas/LESdelta/LESdelta.C
./LES/LESdeltas/cubeRootVolDelta/cubeRootVolDelta.C
./LES/LESdeltas/PrandtlDelta/PrandtlDelta.C
./LES/LESdeltas/vanDriestDelta/vanDriestDelta.C
./LES/LESdeltas/smoothDelta/smoothDelta.C
./LES/LESdeltas/maxDeltaxyz/maxDeltaxyz.C
./LES/LESdeltas/IDDESDelta/IDDESDelta.C
./LES/LESdeltas/IDDESDelta2/IDDESDelta2.C
./LES/LESdeltas/maxDeltaxyzCubeRootLESDelta/maxDeltaxyzCubeRootLESDelta.C
./LES/LESdeltas/maxFaceDelta/maxFaceDelta.C
./LES/LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.C
./LES/LESdeltas/SLADelta/SLADelta.C
./LES/LESfilters/LESfilter/LESfilter.C
./LES/LESfilters/simpleFilter/simpleFilter.C
./LES/LESfilters/laplaceFilter/laplaceFilter.C
./LES/LESfilters/anisotropicFilter/anisotropicFilter.C
./DES/DESModel/DESModelBase.C
./derivedFvPatchFields/fixedShearStress/fixedShearStressFvPatchVectorField.C
./derivedFvPatchFields/porousPressure/porousPressureFvPatchScalarField.C
./derivedFvPatchFields/porousBafflePressure/porousBafflePressureFvPatchField.C
./derivedFvPatchFields/porousBafflePressureAMI/porousBafflePressureAMIFvPatchField.C
./derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.C
./derivedFvPatchFields/relativeHumidityInletOutlet/relativeHumidityInletOutletFvPatchScalarField.C
./derivedFvPatchFields/turbulentMassFluxConcentration/turbulentMassFluxConcentrationFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/analyticalProfiles/FkEBLagAnalyticalProfile/FkEBLagAnalyticalProfile.C
./derivedFvPatchFields/wallFunctions/analyticalProfiles/ePlusEBLagAnalyticalProfile/ePlusEBLagAnalyticalProfile.C
./derivedFvPatchFields/wallFunctions/analyticalProfiles/v2kEBLagAnalyticalProfile/v2kEBLagAnalyticalProfile.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/roughWallCoefficients/roughWallCoefficients.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkWallFunction/nutkWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkRoughWallFunction/nutkRoughWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkGeneralRoughWallFunction/nutkGeneralRoughWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkAtmRoughWallFunction/nutkAtmRoughWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkDemirelRoughWallFunction/nutkDemirelRoughWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutUWallFunction/nutUWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutUSpaldingWallFunction/nutUSpaldingWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutUTabulatedWallFunction/nutUTabulatedWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutLowReWallFunction/nutLowReWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutURoughWallFunction/nutURoughWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutUKnopKOmegaWallFunction/nutUKnopKOmegaWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutUKnopSAWallFunction/nutUKnopSAWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/nutWallFunctions/nutHybridWallFunction/nutHybridWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonEBLagWallFunction/epsilonEBLagWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/omegaWallFunctions/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/omegaWallFunctions/hybridOmegaWallFunction/hybridOmegaWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/omegaWallFunctions/knopOmegaWallFunction/knopOmegaWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/kqRWallFunctions/kqRWallFunction/kqRWallFunctionFvPatchFields.C
./derivedFvPatchFields/wallFunctions/kqRWallFunctions/BlockKqRWallFunction/BlockKqRWallFunctionFvPatchFields.C
./derivedFvPatchFields/wallFunctions/kqRWallFunctions/kLowReWallFunction/kLowReWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/kqRWallFunctions/kEBWallFunction/kEBWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/kqRWallFunctions/kEBLagWallFunction/kEBLagWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/v2WallFunctions/v2WallFunction/v2WallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/v2kWallFunctions/v2kWallFunction/v2kWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/fWallFunctions/fWallFunction/fWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/aWallFunctions/aWallFunction/aWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/atmOmegaWallFunction/atmOmegaWallFunctionFvPatchScalarField.C
./derivedFvPatchFields/wallFunctions/atmEpsilonWallFunction/atmEpsilonWallFunctionFvPatchScalarField.C
./RAS/derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
./RAS/derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
./RAS/derivedFvPatchFields/atmBoundaryLayer/atmBoundaryLayer.C
./RAS/derivedFvPatchFields/atmBoundaryLayerInletVelocity/atmBoundaryLayerInletVelocityFvPatchVectorField.C
./RAS/derivedFvPatchFields/atmBoundaryLayerInletK/atmBoundaryLayerInletKFvPatchScalarField.C
./RAS/derivedFvPatchFields/atmBoundaryLayerInletEpsilon/atmBoundaryLayerInletEpsilonFvPatchScalarField.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        incompressibleTurbulenceModels-includes
    PUBLIC
        fluidThermophysicalModels-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE materialModels
    PRIVATE finiteVolume
    PRIVATE meshTools
    PUBLIC fluidThermophysicalModels
)
