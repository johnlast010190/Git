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
    (c) 2019-2023 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "fvMotionSolvers")

add_helyx_library(${TARGET_NAME} SHARED
./fvMotionSolvers/fvMotionSolver/fvMotionSolver.C
./fvMotionSolvers/displacement/SBRStress/displacementSBRStressFvMotionSolver.C
./fvMotionSolvers/displacement/laplacian/displacementLaplacianFvMotionSolver.C
./fvMotionSolvers/displacement/surfaceAlignedSBRStress/surfaceAlignedSBRStressFvMotionSolver.C
./fvMotionSolvers/displacement/solidBodyDisplacementLaplacian/solidBodyDisplacementLaplacianFvMotionSolver.C
./fvMotionSolvers/displacement/deformingBodyMotion/deformingBodyFvMotionSolver.C
./fvMotionSolvers/componentDisplacement/componentLaplacian/displacementComponentLaplacianFvMotionSolver.C
./fvMotionSolvers/velocity/laplacian/velocityLaplacianFvMotionSolver.C
./fvMotionSolvers/componentVelocity/componentLaplacian/velocityComponentLaplacianFvMotionSolver.C
./fvMotionSolvers/modes/modesFvMotionSolver.C
./fvMotionSolvers/sinus/sinusFvMotionSolver.C
./motionDiffusivity/motionDiffusivity/motionDiffusivity.C
./motionDiffusivity/uniform/uniformDiffusivity.C
./motionDiffusivity/inverseDistance/inverseDistanceDiffusivity.C
./motionDiffusivity/inverseFaceDistance/inverseFaceDistanceDiffusivity.C
./motionDiffusivity/inversePointDistance/inversePointDistanceDiffusivity.C
./motionDiffusivity/inverseVolume/inverseVolumeDiffusivity.C
./motionDiffusivity/directional/directionalDiffusivity.C
./motionDiffusivity/motionDirectional/motionDirectionalDiffusivity.C
./motionDiffusivity/file/fileDiffusivity.C
./motionDiffusivity/fixedRamp/fixedRampDiffusivity.C
./motionDiffusivity/inverseDistAndVol/inverseDistAndVolDiffusivity.C
./motionDiffusivity/manipulators/quadratic/quadraticDiffusivity.C
./motionDiffusivity/manipulators/exponential/exponentialDiffusivity.C
./motionInterpolation/motionInterpolation/motionInterpolation.C
./motionInterpolation/patchCorrected/patchCorrectedInterpolation.C
./motionInterpolation/patchTransformed/patchTransformedInterpolation.C
./fvPatchFields/derived/cellMotion/cellMotionFvPatchFields.C
./fvPatchFields/derived/surfaceSlipDisplacement/surfaceSlipDisplacementFvPatchFields.C
./pointPatchFields/derived/oscillatingVelocity/oscillatingVelocityPointPatchVectorField.C
./pointPatchFields/derived/sedimentation/sedimentationPointPatchVectorField.C
./pointPatchFields/derived/angularOscillatingVelocity/angularOscillatingVelocityPointPatchVectorField.C
./pointPatchFields/derived/constrainedPressureVelocity/constrainedPressureVelocityPointPatchVectorField.C
./pointPatchFields/derived/oscillatingDisplacement/oscillatingDisplacementPointPatchVectorField.C
./pointPatchFields/derived/constrainedOscillatingDisplacement/constrainedOscillatingDisplacementPointPatchVectorField.C
./pointPatchFields/derived/angularOscillatingDisplacement/angularOscillatingDisplacementPointPatchVectorField.C
./pointPatchFields/derived/surfaceSlipDisplacement/surfaceSlipDisplacementPointPatchVectorField.C
./pointPatchFields/derived/surfaceDisplacement/surfaceDisplacementPointPatchVectorField.C
./pointPatchFields/derived/waveDisplacement/waveDisplacementPointPatchVectorField.C
./pointPatchFields/derived/timeVaryingMappedFixedValue/timeVaryingMappedFixedValuePointPatchFields.C
./pointPatchFields/derived/uniformInterpolatedDisplacement/uniformInterpolatedDisplacementPointPatchVectorField.C
./pointPatchFields/derived/bufferZoneFixedDisplacement/bufferZoneFixedDisplacementPointPatchVectorField.C
./pointPatchFields/derived/bufferZoneFixedVelocity/bufferZoneFixedVelocityPointPatchVectorField.C
./pointPatchFields/derived/bufferZoneDisplacementBase/bufferZoneDisplacementBasePointPatchVectorField.C
./pointPatchFields/derived/modalForcedMotion/modalForcedMotionPointPatchVectorField.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fileFormats-includes
        surfMesh-includes
        forces-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
    PRIVATE dynamicMesh
    PRIVATE dynamicFvMesh
    PRIVATE finiteVolume
    PRIVATE fvMeshMovers
)
