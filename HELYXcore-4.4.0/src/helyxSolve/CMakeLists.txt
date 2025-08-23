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
    (c) 2020-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "helyxSolve")

add_helyx_library(${TARGET_NAME} SHARED

# helyxSolve main controller

./helyxSolve.C

# Solver Objects

./solverObjects/solverObject/solverObject.C
./solverObjects/solverObject/solverObjectList.C
./solverObjects/fvMesh/fvMeshSolver.C
./solverObjects/none/noneSolver.C
./solverObjects/flowSolver/flowSolver.C
./solverObjects/rhieChowFlowSolver/rhieChowFlowSolver.C
./solverObjects/singlePhaseSegregatedFlowSolver/singlePhaseSegregatedFlowSolver.C
./solverObjects/fluidEnergySolver/fluidEnergySolver.C
./solverObjects/fixedTemperature/fixedTemperature.C
./solverObjects/turbulence/turbulenceSolver.C
./solverObjects/turbulence/phaseTurbulenceStabilisation/phaseTurbulenceStabilisation.C
./solverObjects/solidHeatConductionSolver/solidHeatConductionSolver.C
./solverObjects/porousSolidThermalEquilibrium/porousSolidThermalEquilibrium.C
./solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.C
./solverObjects/electrostaticPotentialSolver/electrostaticPotentialSolver.C
./solverObjects/electrostaticPotentialSolver/fvPatchFields/electricalBoundaryBase/electricalBoundaryBase.C
./solverObjects/electrostaticPotentialSolver/fvPatchFields/fixedElectricCurrent/fixedElectricCurrentFvPatchScalarField.C
./solverObjects/electrostaticPotentialSolver/jouleHeatingSource/jouleHeatingSource.C
./solverObjects/radiation/radiationSolver.C
./solverObjects/reactions/reactionsSource.C
./solverObjects/reactions/pyrolysisSource.C
./solverObjects/lagrangianSolver/lagrangianSolver.C
./solverObjects/surfaceFilmSolver/surfaceFilmSolver.C
./solverObjects/driftSedimentationSolver/driftSedimentationSolver.C
./solverObjects/windDrivenRain/windDrivenRain.C
./solverObjects/concentrationTransport/concentrationTransport.C
./solverObjects/concentrationTransport/humidityThermalSource/humidityThermalSource.C
./solverObjects/concentrationTransport/mistThermalAndMassSource/mistThermalAndMassSource.C
./solverObjects/nonParticipatingRadiation/nonParticipatingRadiation.C
./solverObjects/multiphaseSegregatedFlowSolver/multiphaseSegregatedFlowSolver.C
./solverObjects/MULESVolumeFractionSolver/MULESVolumeFractionSolver.C
./solverObjects/multiphaseEulerSegregatedFlowSolver/multiphaseEulerSegregatedFlowSolver.C
./solverObjects/phaseChangeSource/phaseChangeSource.C

./solverObjects/disperseEulerian/disperseEulerian.C
./solverObjects/disperseEulerian/disperseEulerianSource/disperseEulerianSource.C
./solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/massTransferModel/massTransferModel.C
./solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/massTransferModel/massTransferModelNew.C
./solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/physisorption/RanzMarshall/RanzMarshallSh.C
./solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/phaseChange/Lee/Lee.C
./solverObjects/disperseEulerian/disperseEulerianSource/massTransfer/noMassTransfer/noMassTransfer.C
./solverObjects/disperseEulerian/disperseEulerianSource/heatTransfer/heatTransferModel/heatTransferModel.C
./solverObjects/disperseEulerian/disperseEulerianSource/heatTransfer/heatTransferModel/heatTransferModelNew.C
./solverObjects/disperseEulerian/disperseEulerianSource/heatTransfer/RanzMarshall/RanzMarshall.C
./solverObjects/disperseEulerian/disperseEulerianSource/heatTransfer/noHeatTransfer/noHeatTransfer.C
./solverObjects/disperseEulerian/phase/phase.C
./solverObjects/disperseEulerian/forceModels/dragModels/dragModel/dragModel.C
./solverObjects/disperseEulerian/forceModels/dragModels/dragModel/dragModelNew.C
./solverObjects/disperseEulerian/forceModels/dragModels/SchillerNaumann/SchillerNaumann.C
./solverObjects/disperseEulerian/forceModels/dragModels/GunnKinzer/GunnKinzer.C
./solverObjects/disperseEulerian/forceModels/liftModels/liftModel/liftModel.C
./solverObjects/disperseEulerian/forceModels/liftModels/liftModel/liftModelNew.C
./solverObjects/disperseEulerian/forceModels/liftModels/constantLiftCoefficient/constantLiftCoefficient.C
./solverObjects/disperseEulerian/forceModels/liftModels/noLift/noLift.C
./solverObjects/disperseEulerian/forceModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.C
./solverObjects/disperseEulerian/forceModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModelNew.C
./solverObjects/disperseEulerian/forceModels/turbulentDispersionModels/Gosman/Gosman.C
./solverObjects/disperseEulerian/forceModels/turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.C
./solverObjects/pcweAcousticSolver/pcweAcousticSolver.C

# Solver functions

./solverFunctionObjects/solverFunctions.C

# Region-coupling

./regionCoupled/fvPatchFields/regionCoupled/regionCoupledFvPatchFields.C
./regionCoupled/fvPatchFields/regionCoupledFlux/regionCoupledEnergyFluxFvPatchScalarField.C
./regionCoupled/fvPatchFields/regionCoupledFlux/regionCoupledFluxFvPatchScalarField.C
./regionCoupled/fvPatchFields/coupledWallEnergy/coupledWallEnergyFvPatchScalarField.C
./regionCoupled/fvPatchFields/specieTransferMassFraction/specieTransferMassFractionFvPatchScalarField.C
./regionCoupled/fvPatchFields/semiPermeableBaffleMassFraction/semiPermeableBaffleMassFractionFvPatchScalarField.C
./regionCoupled/fvPatchFields/specieTransferVelocity/specieTransferVelocityFvPatchVectorField.C
./regionCoupled/fvPatchFields/specieTransferTemperature/specieTransferTemperatureFvPatchScalarField.C
./regionCoupled/compactGaussLaplacianScheme/compactGaussLaplacianSchemes.C

# Monolithic solve

./monolithicSolve/monolithicSolve.C
./monolithicSolve/offsetLduInterface/offsetLduInterfaceBase.C
./monolithicSolve/offsetLduInterface/offsetLduInterface.C
./monolithicSolve/offsetLduInterface/offsetGAMGInterface.C
./monolithicSolve/offsetLduInterface/offsetLduInterfaceFieldBase.C
./monolithicSolve/offsetLduInterface/offsetLduInterfaceField.C
./monolithicSolve/offsetLduInterface/offsetGAMGInterfaceField.C
./monolithicSolve/lduRegisteredPrimitiveMesh/lduRegisteredPrimitiveMesh.C

# Motion solvers

./fvMotionSolvers/displacement/laplacianPeriodic/displacementLaplacianPeriodicFvMotionSolver.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
    PRIVATE dynamicMesh
    PRIVATE finiteVolume
    PRIVATE fvOptions
    PRIVATE turbulenceModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE compressibleTurbulenceModels
    PUBLIC regionModels
    PRIVATE sampling
    PRIVATE dynamicFvMesh
    PRIVATE topoChangerFvMesh
    PRIVATE materialModels
    PRIVATE radiationModels
    PRIVATE ODE
    PRIVATE combustionModels
    PRIVATE chemistryModel
    PUBLIC turbulenceModelSchemes
    PUBLIC solidThermo
    PUBLIC incompressibleTransportModels
    PUBLIC fluidThermophysicalModels
    PUBLIC multiphaseMaterials
    PUBLIC eulerianMultiphaseSystems
    PUBLIC specie
    PRIVATE fvMotionSolvers
    PRIVATE SLGThermo
    PRIVATE lagrangianIntermediate
    PRIVATE lagrangianTurbulence
    PRIVATE thermophysicalProperties
    PRIVATE lagrangian
    PRIVATE coalCombustion
    PRIVATE lagrangianSpray
)

if(${PRECICE_FOUND})
    add_subdirectory(solverObjects/preciceAdapter)
else()
    add_dummy_target(preciceAdapterSolverObject preCICE)
endif()
