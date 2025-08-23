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

set(TARGET_NAME "regionModels")

add_helyx_library(${TARGET_NAME} SHARED
./regionModel/regionModel/regionModel.C
./regionModel/singleLayerRegion/singleLayerRegion.C
./regionModel/regionProperties/regionProperties.C
./regionModel/solutionRegionMapping/solutionRegionMapping.C
./regionModel/solutionInstance/solutionInstance.C
./regionModel/solutionRegions/solutionScheduler.C

./pyrolysisModels/pyrolysisModel/pyrolysisModel.C
./pyrolysisModels/pyrolysisModel/pyrolysisModelNew.C
./pyrolysisModels/reactingOneDim/reactingOneDim.C
./pyrolysisModels/thermo/pyrolysisThermo.C
./pyrolysisModels/pyrolysisModel/pyrolysisModelCollection.C

./surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.C
./surfaceFilmModels/surfaceFilmModel/surfaceFilmModelNew.C
./surfaceFilmModels/noFilm/noFilm.C
./surfaceFilmModels/kinematicSingleLayer/kinematicSingleLayer.C
./surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.C
./surfaceFilmModels/submodels/filmSubModelBase.C
./surfaceFilmModels/submodels/kinematic/force/force/force.C
./surfaceFilmModels/submodels/kinematic/force/force/forceNew.C
./surfaceFilmModels/submodels/kinematic/force/forceList/forceList.C
./surfaceFilmModels/submodels/kinematic/force/contactAngleForces/contactAngleForce/contactAngleForce.C
./surfaceFilmModels/submodels/kinematic/force/contactAngleForces/distribution/distributionContactAngleForce.C
./surfaceFilmModels/submodels/kinematic/force/contactAngleForces/temperatureDependent/temperatureDependentContactAngleForce.C
./surfaceFilmModels/submodels/kinematic/force/contactAngleForces/perturbedTemperatureDependent/perturbedTemperatureDependentContactAngleForce.C
./surfaceFilmModels/submodels/kinematic/force/thermocapillaryForce/thermocapillaryForce.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModel/ejectionModel.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModel/ejectionModelNew.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/ejectionModelList/ejectionModelList.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/drippingEjection/drippingEjection.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/BrunDrippingEjection/BrunDrippingEjection.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/patchEjection/patchEjection.C
./surfaceFilmModels/submodels/kinematic/ejectionModel/curvatureSeparation/curvatureSeparation.C
./surfaceFilmModels/submodels/kinematic/transferModels/transferModel/transferModel.C
./surfaceFilmModels/submodels/kinematic/transferModels/transferModel/transferModelNew.C
./surfaceFilmModels/submodels/kinematic/transferModels/transferModelList/transferModelList.C
./surfaceFilmModels/surfaceFilmRegionModel/surfaceFilmRegionModel.C
./surfaceFilmModels/submodels/kinematic/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModel.C
./surfaceFilmModels/submodels/kinematic/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModelNew.C
./surfaceFilmModels/submodels/kinematic/filmTurbulenceModel/laminar/laminar.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModelNew.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/noPhaseChange/noPhaseChange.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/speciePhaseChange/speciePhaseChange.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/standardPhaseChange/standardPhaseChange.C
./surfaceFilmModels/submodels/thermo/phaseChangeModel/solidification/solidification.C
./surfaceFilmModels/submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModel.C
./surfaceFilmModels/submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModelNew.C
./surfaceFilmModels/submodels/thermo/heatTransferModel/constantHeatTransfer/constantHeatTransfer.C
./surfaceFilmModels/submodels/thermo/heatTransferModel/mappedConvectiveHeatTransfer/mappedConvectiveHeatTransfer.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModelNew.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/noRadiation/noRadiation.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/constantRadiation/constantRadiation.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/primaryRadiation/primaryRadiation.C
./surfaceFilmModels/submodels/thermo/filmRadiationModel/standardRadiation/standardRadiation.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/filmViscosityModel/filmViscosityModel.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/filmViscosityModel/filmViscosityModelNew.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/constantViscosity/constantViscosity.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/NewtonianViscosity/NewtonianViscosity.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/thixotropicViscosity/thixotropicViscosity.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/ArrheniusViscosity/ArrheniusViscosity.C
./surfaceFilmModels/submodels/kinematic/filmViscosityModel/function1Viscosity/function1Viscosity.C
./surfaceFilmModels/derivedFvPatchFields/filmHeightInletVelocity/filmHeightInletVelocityFvPatchVectorField.C
./surfaceFilmModels/derivedFvPatchFields/inclinedFilmNusseltHeight/inclinedFilmNusseltHeightFvPatchScalarField.C
./surfaceFilmModels/derivedFvPatchFields/inclinedFilmNusseltInletVelocity/inclinedFilmNusseltInletVelocityFvPatchVectorField.C
./surfaceFilmModels/derivedFvPatchFields/wallFunctions/alphatFilmWallFunction/alphatFilmWallFunctionFvPatchScalarField.C
./surfaceFilmModels/derivedFvPatchFields/wallFunctions/nutkFilmWallFunction/nutkFilmWallFunctionFvPatchScalarField.C

./regionCoupling/derivedFvPatchFields/filmPyrolysisVelocityCoupled/filmPyrolysisVelocityCoupledFvPatchVectorField.C
./regionCoupling/derivedFvPatchFields/filmPyrolysisTemperatureCoupled/filmPyrolysisTemperatureCoupledFvPatchScalarField.C
./regionCoupling/derivedFvPatchFields/filmPyrolysisRadiativeCoupledMixed/filmPyrolysisRadiativeCoupledMixedFvPatchScalarField.C
./regionCoupling/derivedFvPatchFields/filmTemperatureCoupledBaffleMixed/filmTemperatureCoupledBaffleMixedFvPatchScalarField.C

./thermalBaffle/thermalBaffle/thermalBaffle.C
./thermalBaffle/derivedFvPatchFields/thermalBaffle/thermalBaffleFvPatchScalarField.C
)


helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        chemistryModel-includes
        specie-includes
        fluidThermophysicalModels-includes
        solidThermo-includes
        thermophysicalProperties-includes
        radiationModels-includes
        specie-includes
        compressibleTurbulenceModels-includes
        extrudeModel-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE dynamicFvMesh
    PRIVATE dynamicMesh
    PRIVATE chemistryModel
    PRIVATE materialModels
    PRIVATE specie
    PRIVATE fluidThermophysicalModels
    PRIVATE solidThermo
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE radiationModels
    PRIVATE thermophysicalProperties
    PRIVATE distributionModels
    PRIVATE multiphaseMaterials
)


