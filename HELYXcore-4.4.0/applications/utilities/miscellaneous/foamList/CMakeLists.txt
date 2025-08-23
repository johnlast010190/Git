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

set(TARGET_NAME "foamList-exe")

add_helyx_executable(${TARGET_NAME}
./foamList.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/transportModels
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE barotropicCompressibilityModel
    PRIVATE blockMesh
    PRIVATE chemistryModel
    PRIVATE coalCombustion
    PRIVATE combustionModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE conversion
    PRIVATE decompose
    PRIVATE decompositionMethods
    PRIVATE distributed
    PRIVATE distributionModels
    PRIVATE DPMTurbulenceModels
    PRIVATE driftFluxRelativeVelocityModels
    PRIVATE driftFluxTransportModels
    PRIVATE DSMC
    PRIVATE dynamicFvMesh
    PRIVATE dynamicMesh
    PRIVATE engine
    PRIVATE extrudeModel
    PRIVATE fieldFunctionObjects
    PRIVATE fileFormats
    PRIVATE finiteVolume
    PRIVATE fluidThermophysicalModels
    PRIVATE forces
    PRIVATE fvMotionSolvers
    PRIVATE fvOptions
    PRIVATE genericPatchFields
    PRIVATE incompressibleTransportModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE lagrangian
    PRIVATE lagrangianFunctionObjects
    PRIVATE lagrangianIntermediate
    PRIVATE lagrangianSpray
    PRIVATE lagrangianTurbulence
    PRIVATE laminarFlameSpeedModels
    PRIVATE meshTools
    PRIVATE molecularMeasurements
    PRIVATE molecule
    PRIVATE ODE
    PRIVATE pairPatchAgglomeration
    PRIVATE potential
    PRIVATE radiationModels
    PRIVATE reconstruct
    PRIVATE regionModels
    PRIVATE renumberMethods
    PRIVATE rhoCentralFoam
    PRIVATE rigidBodyDynamics
    PRIVATE rigidBodyMeshMotion
    PRIVATE sampling
    PRIVATE ${SCOTCHDECOMP}
    PRIVATE sixDoFRigidBodyMotion
    PRIVATE SLGThermo
    PRIVATE snappyHexMesh
    PRIVATE solidParticle
    PRIVATE solidThermo
    PRIVATE helyxSolve
    PRIVATE specie
    PRIVATE surfMesh
    PRIVATE thermophysicalProperties
    PRIVATE topoChangerFvMesh
    PRIVATE turbulenceModels
    PRIVATE twoPhaseMixture
    PRIVATE utilityFunctionObjects
    PRIVATE helyxSolve
)
