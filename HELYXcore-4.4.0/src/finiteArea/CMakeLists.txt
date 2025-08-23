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

set(TARGET_NAME "finiteArea")

add_helyx_library(${TARGET_NAME} SHARED
./faMesh/faGlobalMeshData/faGlobalMeshData.C
./faMesh/faMesh.C
./faMesh/faMeshDemandDrivenData.C
./faMesh/faMeshUpdate.C
./faMesh/faBoundaryMesh/faBoundaryMesh.C
./faMesh/faPatches/faPatch/faPatch.C
./faMesh/faPatches/faPatch/newFaPatch.C
./faMesh/faPatches/basic/coupled/coupledFaPatch.C
./faMesh/faPatches/constraint/empty/emptyFaPatch.C
./faMesh/faPatches/constraint/processor/processorFaPatch.C
./faMesh/faPatches/constraint/wedge/wedgeFaPatch.C
./faMesh/faPatches/constraint/cyclic/cyclicFaPatch.C
./faMesh/faPatches/constraint/symmetry/symmetryFaPatch.C
./faMesh/faMeshMapper/faMeshMapper.C
./faMesh/faMeshMapper/faAreaMapper.C
./faMesh/faMeshMapper/faEdgeMapper.C
./faMesh/faMeshMapper/faPatchMapper.C
./fields/faPatchFields/faPatchField/faPatchFields.C
./fields/faPatchFields/basic/basicSymmetry/basicSymmetryFaPatchScalarField.C
./fields/faPatchFields/basic/calculated/calculatedFaPatchFields.C
./fields/faPatchFields/basic/coupled/coupledFaPatchFields.C
./fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.C
./fields/faPatchFields/basic/fixedValue/fixedValueFaPatchFields.C
./fields/faPatchFields/basic/fixedGradient/fixedGradientFaPatchFields.C
./fields/faPatchFields/basic/mixed/mixedFaPatchFields.C
./fields/faPatchFields/basic/transform/transformFaPatchFields.C
./fields/faPatchFields/basic/transform/transformFaPatchScalarField.C
./fields/faPatchFields/constraint/empty/emptyFaPatchFields.C
./fields/faPatchFields/constraint/processor/processorFaPatchFields.C
./fields/faPatchFields/constraint/processor/processorFaPatchScalarField.C
./fields/faPatchFields/constraint/wedge/wedgeFaPatchFields.C
./fields/faPatchFields/constraint/wedge/wedgeFaPatchScalarField.C
./fields/faPatchFields/constraint/cyclic/cyclicFaPatchFields.C
./fields/faPatchFields/constraint/symmetry/symmetryFaPatchFields.C
./fields/faPatchFields/derived/fixedValueOutflow/fixedValueOutflowFaPatchFields.C
./fields/faPatchFields/derived/inletOutlet/inletOutletFaPatchFields.C
./fields/faPatchFields/derived/slip/slipFaPatchFields.C
./fields/faPatchFields/derived/edgeNormalFixedValue/edgeNormalFixedValueFaPatchVectorField.C
./fields/faPatchFields/derived/uniformFixedValue/uniformFixedValueFaPatchFields.C
./fields/faPatchFields/derived/fixedValueZeroGradient/fixedValueZeroGradientFaPatchFields.C
./faMesh/geodeticWallDist/geodeticWallDist/geodeticWallDist.C
./faMesh/geodeticWallDist/faPatchDistMethods/faPatchDistMethod/faPatchDistMethod.C
./faMesh/geodeticWallDist/faPatchDistMethods/Poisson/PoissonFaPatchDistMethod.C
./faMesh/geodeticWallDist/faPatchDistMethods/advectionDiffusion/advectionDiffusionFaPatchDistMethod.C
./fields/faePatchFields/faePatchField/faePatchFields.C
./fields/faePatchFields/basic/calculated/calculatedFaePatchFields.C
./fields/faePatchFields/basic/coupled/coupledFaePatchFields.C
./fields/faePatchFields/basic/fixedValue/fixedValueFaePatchFields.C
./fields/faePatchFields/constraint/empty/emptyFaePatchFields.C
./fields/faePatchFields/constraint/processor/processorFaePatchFields.C
./fields/faePatchFields/constraint/wedge/wedgeFaePatchFields.C
./fields/faePatchFields/constraint/cyclic/cyclicFaePatchFields.C
./fields/faePatchFields/constraint/symmetry/symmetryFaePatchFields.C
./fields/areaFields/areaFields.C
./fields/edgeFields/edgeFields.C
./faMatrices/faMatrices.C
./faMatrices/faScalarMatrix/faScalarMatrix.C
./interpolation/edgeInterpolation/edgeInterpolation.C
./interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationSchemes.C
./interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolationMake.C
./interpolation/edgeInterpolation/schemes/upwind/upwindEdgeInterpolationMake.C
./interpolation/edgeInterpolation/schemes/linearUpwind/linearUpwindEdgeInterpolationMake.C
./interpolation/edgeInterpolation/schemes/Gamma/GammaEdgeInterpolationMake.C
./interpolation/edgeInterpolation/schemes/blended/blendedEdgeInterpolationMake.C
./finiteArea/fa/fa.C
./finiteArea/faSchemes/faSchemes.C
./finiteArea/ddtSchemes/faDdtScheme/faDdtSchemes.C
./finiteArea/ddtSchemes/steadyStateFaDdtScheme/steadyStateFaDdtSchemes.C
./finiteArea/ddtSchemes/EulerFaDdtScheme/EulerFaDdtSchemes.C
./finiteArea/ddtSchemes/backwardFaDdtScheme/backwardFaDdtSchemes.C
./finiteArea/ddtSchemes/boundedBackwardFaDdtScheme/boundedBackwardFaDdtScheme.C
./finiteArea/d2dt2Schemes/faD2dt2Scheme/faD2dt2Schemes.C
./finiteArea/d2dt2Schemes/EulerFaD2dt2Scheme/EulerFaD2dt2Schemes.C
./finiteArea/divSchemes/faDivScheme/faDivSchemes.C
./finiteArea/divSchemes/gaussFaDivScheme/gaussFaDivSchemes.C
./finiteArea/gradSchemes/faGradScheme/faGradSchemes.C
./finiteArea/gradSchemes/gaussFaGrad/gaussFaGrads.C
./finiteArea/gradSchemes/leastSquaresFaGrad/leastSquaresFaVectors.C
./finiteArea/gradSchemes/leastSquaresFaGrad/leastSquaresFaGrads.C
./finiteArea/gradSchemes/limitedGradSchemes/faceLimitedFaGrad/faceLimitedFaGrads.C
./finiteArea/gradSchemes/limitedGradSchemes/edgeLimitedFaGrad/edgeLimitedFaGrads.C
./finiteArea/lnGradSchemes/lnGradScheme/lnGradSchemes.C
./finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrads.C
./finiteArea/lnGradSchemes/limitedLnGrad/limitedLnGrads.C
./finiteArea/lnGradSchemes/fourthLnGrad/fourthLnGrads.C
./finiteArea/lnGradSchemes/uncorrectedLnGrad/uncorrectedLnGrads.C
./finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianSchemes.C
./finiteArea/laplacianSchemes/gaussFaLaplacianScheme/gaussFaLaplacianSchemes.C
./finiteArea/convectionSchemes/faConvectionScheme/faConvectionSchemes.C
./finiteArea/convectionSchemes/gaussFaConvectionScheme/gaussFaConvectionSchemes.C
./finiteArea/fac/facMeshPhi.C
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/interpolation/edgeInterpolation/schemes/Gamma/GammaEdgeInterpolationMake.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        finiteVolume-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE OpenFOAM
    PRIVATE meshTools
)
