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
    (c) 2019-2024 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "OpenFOAM")

# Set the current git hash into global.C
# We want this to run whenever the current git hash is changed, so we use
# the timestamp of .git/index to detect any change in the git state
if(EXISTS "${HELYX_PROJECT_DIR}/.git/index")
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/global.C
        COMMAND
            ${CMAKE_COMMAND} -DinputFile=./global/global.C.template -DoutputFile=${CMAKE_CURRENT_BINARY_DIR}/global.C -DHELYX_PROJECT_VERSION=${HELYX_PROJECT_VERSION} -DHELYX_PROJECT_DIR=${HELYX_PROJECT_DIR} -P ${CMAKE_SOURCE_DIR}/etc/cmake/scripts/buildTimeSetGitHash.cmake
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/global/global.C.template
        DEPENDS ${HELYX_PROJECT_DIR}/.git/index
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
else()
    message(STATUS ".git directory not found in \$HELYX_PROJECT_DIR, unable to find git hash")
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/global.C
        COMMAND
            ${CMAKE_COMMAND} -DinputFile=./global/global.C.template -DoutputFile=${CMAKE_CURRENT_BINARY_DIR}/global.C -DGIT_HASH="unknown" -DHELYX_PROJECT_VERSION=${HELYX_PROJECT_VERSION} -DHELYX_PROJECT_DIR=${HELYX_PROJECT_DIR} -P ${CMAKE_SOURCE_DIR}/etc/cmake/scripts/buildTimeFileConfigurer.cmake
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/global/global.C.template
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# Fix issues with FPE for tree vectorisation run in SP on AMD
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "$ENV{HELYX_PRECISION_OPTION}" MATCHES "SP")
file(GLOB NOTREE_SRC_FILES "./meshes/meshShapes/cell/cell.C")
set_property(SOURCE ${NOTREE_SRC_FILES} PROPERTY COMPILE_FLAGS -fno-tree-vectorize)
endif()

# Don't use add_helyx_library here, because this and ptscotchDecomp have non-
# standard output directories.  Because it's only two libraries, it's probably
# not worth developing a general method.
add_library(${TARGET_NAME} SHARED
${CMAKE_CURRENT_BINARY_DIR}/global.C
./global/argList/argList.C
./global/clock/clock.C
./global/debug/simpleObjectRegistry.C
./global/columnatedFileWriter/columnatedFileWriter.C
./global/etcFiles/etcFiles.C
./global/fileOperations/fileOperation/fileOperation.C
./global/fileOperations/fileOperationInitialise/fileOperationInitialise.C
./global/fileOperations/uncollatedFileOperation/uncollatedFileOperation.C
./global/fileOperations/masterUncollatedFileOperation/masterUncollatedFileOperation.C
./global/fileOperations/collatedFileOperation/collatedFileOperation.C
./global/fileOperations/collatedFileOperation/hostCollatedFileOperation.C
./global/fileOperations/collatedFileOperation/threadedCollatedOFstream.C
./global/fileOperations/collatedFileOperation/OFstreamCollator.C
./global/libmFunctions/libmFunctions.C
./primitives/bools/bool/bool.C
./primitives/bools/bool/boolIO.C
./primitives/bools/Switch/Switch.C
./primitives/bools/Switch/SwitchIO.C
./primitives/bools/lists/boolList.C
./primitives/chars/char/charIO.C
./primitives/chars/wchar/wcharIO.C
./primitives/direction/directionIO.C
./primitives/ints/uint32/uint32.C
./primitives/ints/uint32/uint32IO.C
./primitives/ints/uint64/uint64.C
./primitives/ints/uint64/uint64IO.C
./primitives/ints/uint/uintIO.C
./primitives/ints/int32/int32.C
./primitives/ints/int32/int32IO.C
./primitives/ints/int64/int64.C
./primitives/ints/int64/int64IO.C
./primitives/ints/int/intIO.C
./primitives/ints/label/label.C
./primitives/ints/uLabel/uLabel.C
./primitives/ints/lists/labelIOList.C
./primitives/ints/lists/labelListIOList.C
./primitives/Scalar/doubleScalar/doubleScalar.C
./primitives/Scalar/floatScalar/floatScalar.C
./primitives/Scalar/scalar/scalar.C
./primitives/Scalar/scalar/invIncGamma.C
./primitives/Scalar/lists/scalarList.C
./primitives/Scalar/lists/scalarIOList.C
./primitives/Scalar/lists/scalarListIOList.C

# For the unity build, it's necessary to compile vectors before tensors
./primitives/Vector/complexVector/complexVector.C
./primitives/Vector/floatVector/floatVector.C
./primitives/Vector/labelVector/labelVector.C
./primitives/Vector/vector/vector.C
./primitives/Vector/lists/vectorList.C
./primitives/Vector/lists/vectorIOList.C
./primitives/Vector/lists/vectorListIOList.C
./primitives/Vector2D/vector2D/vector2D.C

./primitives/DiagTensor/diagTensor/diagTensor.C
./primitives/SphericalTensor/sphericalTensor/sphericalTensor.C
./primitives/SphericalTensor/labelSphericalTensor/labelSphericalTensor.C
./primitives/SymmTensor/labelSymmTensor/labelSymmTensor.C
./primitives/SymmTensor/symmTensor/symmTensor.C
./primitives/Tensor/labelTensor/labelTensor.C
./primitives/Tensor/tensor/tensor.C
./primitives/Tensor/lists/sphericalTensorList.C
./primitives/Tensor/lists/symmTensorList.C
./primitives/Tensor/lists/tensorList.C
./primitives/Tensor/modifiedNewtonsMethod/modifiedNewtonsMethod.C
./primitives/Tensor4thRank/tensor4thRank/tensor4thRank.C
./primitives/Tensor/floatTensor/floatTensor.C
./primitives/Tensor2D/tensor2D/tensor2D.C
./primitives/SphericalTensor2D/sphericalTensor2D/sphericalTensor2D.C
./primitives/SymmTensor2D/symmTensor2D/symmTensor2D.C
./primitives/complex/complex.C
./primitives/globalIndexAndTransform/globalIndexAndTransform.C
./primitives/transform/transformer/transformer.C
./primitives/quaternion/quaternion.C
./primitives/septernion/septernion.C
./primitives/triad/triad.C
./primitives/functions/Function1/makeDataEntries.C
./primitives/functions/Function1/ramp/ramp.C
./primitives/functions/Function1/linearRamp/linearRamp.C
./primitives/functions/Function1/quadraticRamp/quadraticRamp.C
./primitives/functions/Function1/quarterSineRamp/quarterSineRamp.C
./primitives/functions/Function1/quarterCosineRamp/quarterCosineRamp.C
./primitives/functions/Function1/halfCosineRamp/halfCosineRamp.C
./primitives/functions/Function1/AIJWindProfile/AIJWindProfile.C
./primitives/functions/Function1/AIJWindProfile/AIJWindKineticEnergyProfile.C
./primitives/functions/Function1/AIJWindProfile/AIJWindDissipationRateProfile.C
./primitives/functions/Function1/AIJWindProfile/AIJWindSpecificDissipationRateProfile.C
./primitives/functions/Function1/RHWindProfile/RHWindProfile.C
./primitives/functions/Function1/RHWindProfile/RHWindKineticEnergyProfile.C
./primitives/functions/Function1/RHWindProfile/RHWindDissipationRateProfile.C
./primitives/functions/Function1/RHWindProfile/RHWindSpecificDissipationRateProfile.C
./primitives/functions/Function1/Constant/ConstantSpecialisations.C
./primitives/functions/Function1/NonUniformTable/NonUniformTable1Specialisations.C
./primitives/functions/Function1/Table/TableSpecialisations.C
./primitives/functions/Function1/powerLaw/powerLaw.C
./primitives/functions/Function1/Andrade/Andrade.C
./primitives/functions/Function2/makeFunction2s.C
./primitives/functions/Function2/APIdiffCoef/APIdiffCoef.C
./primitives/functions/Polynomial/polynomialFunction.C
./primitives/subModelBase/subModelBase.C
./primitives/nil/nil.C
./primitives/strings/string/string.C
./primitives/strings/string/stringIO.C
./primitives/strings/string/stringIOList.C
./primitives/strings/word/word.C
./primitives/strings/word/wordIO.C
./primitives/strings/word/wordIOList.C
./primitives/strings/fileName/fileName.C
./primitives/strings/fileName/fileNameIO.C
./primitives/strings/keyType/keyType.C
./primitives/strings/wordRe/wordRe.C
./primitives/strings/wordRes/wordRes.C
./primitives/strings/lists/hashedWordList.C
./primitives/strings/stringOps/stringOps.C
./primitives/strings/stringOps/stringOpsSort.C
./primitives/ops/flipOp.C
./primitives/hashes/Hasher/Hasher.C
./primitives/hashes/SHA1/SHA1.C
./primitives/hashes/SHA1/SHA1Digest.C
./primitives/random/Random/Random.C
./primitives/ranges/labelRange/labelRange.C
./primitives/ranges/labelRange/labelRanges.C
./primitives/ranges/scalarRange/scalarRange.C
./primitives/ranges/scalarRange/scalarRanges.C
./primitives/ranges/sliceRange/sliceRange.C
./primitives/ranges/tableBounds/tableBounds.C
./primitives/spatialVectorAlgebra/SpatialVector/spatialVector/spatialVector.C
./primitives/spatialVectorAlgebra/SpatialTensor/spatialTensor/spatialTensor.C
./primitives/spatialVectorAlgebra/CompactSpatialTensor/compactSpatialTensor/compactSpatialTensor.C
./primitives/polynomialEqns/cubicEqn/cubicEqn.C
./primitives/polynomialEqns/quadraticEqn/quadraticEqn.C
./primitives/Barycentric2D/barycentric2D/barycentric2D.C
./containers/HashTables/HashTable/HashTableCore.C
./containers/HashTables/StaticHashTable/StaticHashTableCore.C
./containers/Lists/SortableList/ParSortableListName.C
./containers/Lists/PackedList/PackedListCore.C
./containers/Lists/PackedList/PackedBoolList.C
./containers/Lists/ListOps/ListOps.C
./containers/LinkedLists/linkTypes/SLListBase/SLListBase.C
./containers/LinkedLists/linkTypes/DLListBase/DLListBase.C
./db/IOstreams/token/tokenIO.C
./db/IOstreams/IOstreams/IOstream.C
./db/IOstreams/IOstreams/Istream.C
./db/IOstreams/IOstreams/Ostream.C
./db/IOstreams/Sstreams/ISstream.C
./db/IOstreams/Sstreams/OSstream.C
./db/IOstreams/Sstreams/SstreamsPrint.C
./db/IOstreams/Sstreams/readHexLabel.C
./db/IOstreams/Sstreams/prefixOSstream.C
./db/IOstreams/hashes/base64Layer.C
./db/IOstreams/gzstream/gzstream.C
./db/IOstreams/Fstreams/IFstream.C
./db/IOstreams/Fstreams/OFstream.C
./db/IOstreams/Fstreams/masterOFstream.C
./db/IOstreams/Tstreams/ITstream.C
./db/IOstreams/StringStreams/StringStream.C
./db/IOstreams/Pstreams/UIPstream.C
./db/IOstreams/Pstreams/IPstream.C
./db/IOstreams/Pstreams/mpiFunctors.C
./db/IOstreams/Pstreams/UPstreamCommsStruct.C
./db/IOstreams/Pstreams/Pstream.C
./db/IOstreams/Pstreams/UOPstream.C
./db/IOstreams/Pstreams/OPstream.C
./db/IOstreams/Pstreams/PstreamBuffers.C
./db/IOstreams/Pstreams/UOPwrite.C
./db/IOstreams/Pstreams/UIPread.C
./db/IOstreams/Pstreams/UPstream.C
./db/IOstreams/Pstreams/PstreamGlobals.C
./db/dictionary/dictionary.C
./db/dictionary/dictionaryIO.C
./db/dictionary/dictionarySearch.C
./db/dictionary/dictionaryCompat.C
./db/dictionary/entry/entry.C
./db/dictionary/entry/entryIO.C
./db/dictionary/primitiveEntry/primitiveEntry.C
./db/dictionary/primitiveEntry/primitiveEntryIO.C
./db/dictionary/dictionaryEntry/dictionaryEntry.C
./db/dictionary/dictionaryEntry/dictionaryEntryIO.C
./db/dictionary/dictionaryListEntry/dictionaryListEntry.C
./db/dictionary/dictionaryListEntry/dictionaryListEntryIO.C
./db/dictionary/functionEntries/calcEntry/calcEntry.C
./db/dictionary/functionEntries/codeStream/codeStream.C
./db/dictionary/functionEntries/functionEntry/functionEntry.C
./db/dictionary/functionEntries/includeEntry/includeEntry.C
./db/dictionary/functionEntries/includeEtcEntry/includeEtcEntry.C
./db/dictionary/functionEntries/includeFuncEntry/includeFuncEntry.C
./db/dictionary/functionEntries/includeIfPresentEntry/includeIfPresentEntry.C
./db/dictionary/functionEntries/inputModeEntry/inputModeEntry.C
./db/dictionary/functionEntries/removeEntry/removeEntry.C
./db/IOobjects/IOdictionary/IOdictionary.C
./db/IOobjects/IOdictionary/IOdictionaryIO.C
./db/IOobjects/IOdictionary/localIOdictionary.C
./db/IOobjects/IOdictionary/systemDict.C
./db/IOobjects/IOMap/IOMapName.C
./db/IOobjects/decomposedBlockData/decomposedBlockData.C
./db/IOobjects/GlobalIOField/GlobalIOFields.C
./db/IOobjects/GlobalIOList/globalIOLists.C
./db/IOobject/IOobject.C
./db/IOobject/IOobjectIO.C
./db/IOobject/IOobjectReadHeader.C
./db/IOobject/IOobjectWriteHeader.C
./db/regIOobject/regIOobjectRead.C
./db/regIOobject/regIOobjectWrite.C
./db/IOobjectList/IOobjectList.C
./db/objectRegistry/objectRegistry.C
./db/CallbackRegistry/CallbackRegistryName.C
./db/dynamicLibrary/dlLibraryTable/dlLibraryTable.C
./db/dynamicLibrary/dynamicCode/dynamicCode.C
./db/dynamicLibrary/dynamicCode/dynamicCodeContext.C
./db/dynamicLibrary/codedBase/codedBase.C
./db/solutionInstanceRegistry/solutionInstanceRegistry.C
./db/solutionRegistry/solutionRegistry.C
./db/functionObjects/functionObject/functionObject.C
./db/functionObjects/nullFunctionObject/nullFunctionObject.C
./db/functionObjects/functionObjectList/functionObjectList.C
./db/functionObjects/stateFunctionObject/stateFunctionObject.C
./db/functionObjects/writeFile/writeFile.C
./db/functionObjects/logFiles/logFiles.C
./db/functionObjects/timeControl/timeControl.C
./db/functionObjects/timeControl/timeControlFunctionObject.C
./db/functionObjects/regionFunctionObject/regionFunctionObject.C
./db/functionObjects/timeFunctionObject/timeFunctionObject.C
./db/Time/TimePaths.C
./db/Time/TimeState.C
./db/Time/Time.C
./db/Time/TimeIO.C
./db/Time/findTimes.C
./db/Time/subCycleTime.C
./db/Time/findInstance.C
./db/Time/timeSelector.C
./db/Time/instant/instant.C
./dimensionSet/dimensionSet.C
./dimensionSet/dimensionSetIO.C
./dimensionSet/dimensionSets.C
./dimensionedTypes/dimensionedScalar/dimensionedScalar.C
./dimensionedTypes/dimensionedSphericalTensor/dimensionedSphericalTensor.C
./dimensionedTypes/dimensionedSymmTensor/dimensionedSymmTensor.C
./dimensionedTypes/dimensionedTensor/dimensionedTensor.C
./orientedType/orientedType.C
./matrices/solution/solution.C
./matrices/scalarMatrices/scalarMatrices.C
./matrices/scalarMatrices/SVD/SVD.C
./matrices/LUscalarMatrix/LUscalarMatrix.C
./matrices/LUscalarMatrix/procLduMatrix.C
./matrices/LUscalarMatrix/procLduInterface.C
./matrices/lduMatrix/lduMatrix/lduMatrix.C
./matrices/lduMatrix/lduMatrix/lduMatrixOperations.C
./matrices/lduMatrix/lduMatrix/lduMatrixATmul.C
./matrices/lduMatrix/lduMatrix/lduMatrixUpdateMatrixInterfaces.C

# For the Unity build, LduMatrix/solverPerformance.C must be included before
# lduMatrix/lduMatrixSolver.C to ensure templates are instatantiated in the
# right order
./matrices/LduMatrix/LduMatrix/solverPerformance.C

./matrices/lduMatrix/lduMatrix/lduMatrixSolver.C
./matrices/lduMatrix/lduMatrix/lduMatrixSmoother.C
./matrices/lduMatrix/lduMatrix/lduMatrixPreconditioner.C
./matrices/lduMatrix/solvers/diagonalSolver/diagonalSolver.C
./matrices/lduMatrix/solvers/smoothSolver/smoothSolver.C
./matrices/lduMatrix/solvers/PCG/PCG.C
./matrices/lduMatrix/solvers/PBiCG/PBiCG.C
./matrices/lduMatrix/solvers/PBiCGStab/PBiCGStab.C
./matrices/lduMatrix/smoothers/GaussSeidel/GaussSeidelSmoother.C
./matrices/lduMatrix/smoothers/symGaussSeidel/symGaussSeidelSmoother.C
./matrices/lduMatrix/smoothers/nonBlockingGaussSeidel/nonBlockingGaussSeidelSmoother.C
./matrices/lduMatrix/smoothers/DIC/DICSmoother.C
./matrices/lduMatrix/smoothers/FDIC/FDICSmoother.C
./matrices/lduMatrix/smoothers/DICGaussSeidel/DICGaussSeidelSmoother.C
./matrices/lduMatrix/smoothers/DILU/DILUSmoother.C
./matrices/lduMatrix/smoothers/DILUGaussSeidel/DILUGaussSeidelSmoother.C
./matrices/lduMatrix/preconditioners/noPreconditioner/noPreconditioner.C
./matrices/lduMatrix/preconditioners/diagonalPreconditioner/diagonalPreconditioner.C
./matrices/lduMatrix/preconditioners/DICPreconditioner/DICPreconditioner.C
./matrices/lduMatrix/preconditioners/FDICPreconditioner/FDICPreconditioner.C
./matrices/lduMatrix/preconditioners/DILUPreconditioner/DILUPreconditioner.C
./matrices/lduMatrix/preconditioners/GAMGPreconditioner/GAMGPreconditioner.C
./matrices/lduMatrix/lduAddressing/lduAddressing.C
./matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.C
./matrices/lduMatrix/lduAddressing/lduInterface/processorLduInterface.C
./matrices/lduMatrix/lduAddressing/lduInterface/cyclicLduInterface.C
./matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.C
./matrices/lduMatrix/lduAddressing/lduInterfaceFields/processorLduInterfaceField/processorLduInterfaceField.C
./matrices/lduMatrix/lduAddressing/lduInterfaceFields/cyclicLduInterfaceField/cyclicLduInterfaceField.C
./matrices/lduMatrix/solvers/GAMG/GAMGSolver.C
./matrices/lduMatrix/solvers/GAMG/GAMGSolverAgglomerateMatrix.C
./matrices/lduMatrix/solvers/GAMG/GAMGSolverInterpolate.C
./matrices/lduMatrix/solvers/GAMG/GAMGSolverScale.C
./matrices/lduMatrix/solvers/GAMG/GAMGSolverSolve.C
./matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.C
./matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterfaceNew.C
./matrices/lduMatrix/solvers/GAMG/interfaces/processorGAMGInterface/processorGAMGInterface.C
./matrices/lduMatrix/solvers/GAMG/interfaces/processorCyclicGAMGInterface/processorCyclicGAMGInterface.C
./matrices/lduMatrix/solvers/GAMG/interfaces/cyclicGAMGInterface/cyclicGAMGInterface.C
./matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.C
./matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceFieldNew.C
./matrices/lduMatrix/solvers/GAMG/interfaceFields/processorGAMGInterfaceField/processorGAMGInterfaceField.C
./matrices/lduMatrix/solvers/GAMG/interfaceFields/processorCyclicGAMGInterfaceField/processorCyclicGAMGInterfaceField.C
./matrices/lduMatrix/solvers/GAMG/interfaceFields/cyclicGAMGInterfaceField/cyclicGAMGInterfaceField.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/GAMGAgglomeration/GAMGAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/GAMGAgglomeration/GAMGAgglomerateLduAddressing.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/pairGAMGAgglomeration/pairGAMGAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/pairGAMGAgglomeration/pairGAMGAgglomerate.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/algebraicPairGAMGAgglomeration/algebraicPairGAMGAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/dummyAgglomeration/dummyAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/GAMGProcAgglomeration/GAMGProcAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/masterCoarsestGAMGProcAgglomeration/masterCoarsestGAMGProcAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/manualGAMGProcAgglomeration/manualGAMGProcAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/eagerGAMGProcAgglomeration/eagerGAMGProcAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/noneGAMGProcAgglomeration/noneGAMGProcAgglomeration.C
./matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/procFacesGAMGProcAgglomeration/procFacesGAMGProcAgglomeration.C
./meshes/lduMesh/lduMesh.C
./meshes/lduMesh/lduPrimitiveMesh.C
./matrices/LduMatrix/LduMatrix/lduMatrices.C
./matrices/LduMatrix/LduMatrix/LduInterfaceField/LduInterfaceFields.C
./matrices/LduMatrix/Smoothers/lduSmoothers.C
./matrices/LduMatrix/Preconditioners/lduPreconditioners.C
./matrices/LduMatrix/Solvers/lduSolvers.C
./memory/tmp/tmp.C
./meshes/primitiveShapes/line/line.C
./meshes/primitiveShapes/plane/plane.C
./meshes/primitiveShapes/triangle/intersection.C
./meshes/primitiveShapes/objectHit/pointIndexHitIOList.C
./meshes/meshShapes/edge/edge.C
./meshes/meshShapes/edge/edgeIOList.C
./meshes/meshShapes/face/face.C
./meshes/meshShapes/face/faceIntersection.C
./meshes/meshShapes/face/faceContactSphere.C
./meshes/meshShapes/face/faceAreaInContact.C
./meshes/meshShapes/face/faceIOList.C
./meshes/meshShapes/cell/cell.C
./meshes/meshShapes/cell/oppositeCellFace.C
./meshes/meshShapes/cell/cellIOList.C
./meshes/meshShapes/tetCell/tetCell.C
./meshes/meshShapes/cellModeller/cellModeller.C
./meshes/meshShapes/cellModel/cellModel.C
./meshes/meshShapes/cellModel/cellModelIO.C
./meshes/meshShapes/cellShape/cellShape.C
./meshes/meshShapes/cellShape/cellShapeEqual.C
./meshes/meshShapes/cellShape/cellShapeIO.C
./meshes/meshShapes/cellShape/cellShapeIOList.C
./meshes/Identifiers/patch/patchIdentifier.C
./meshes/Identifiers/patch/coupleGroupIdentifier.C
./meshes/Identifiers/surface/surfZoneIdentifier.C
./meshes/MeshObject/meshObject.C

# For the Unity build, pointPatchFields.C must be included before polyPatch.C
# to ensure templates are instatantiated in the right order
./fields/pointPatchFields/pointPatchField/pointPatchFields.C

./meshes/polyMesh/polyPatches/polyPatch/polyPatch.C
./meshes/polyMesh/polyPatches/polyPatch/polyPatchNew.C
./meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.C
./meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.C
./meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.C
./meshes/polyMesh/polyPatches/basic/generic/genericPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/cyclic/cyclicTransform.C
./meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/cyclicSlip/cyclicSlipPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/processorCyclic/processorCyclicPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/symmetryPlane/symmetryPlanePolyPatch.C
./meshes/polyMesh/polyPatches/constraint/symmetry/symmetryPolyPatch.C
./meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C
./meshes/polyMesh/polyPatches/constraint/internal/internalPolyPatch.C
./meshes/polyMesh/polyPatches/derived/baseWall/wall.C
./meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.C
./meshes/polyMesh/polyPatches/derived/indirectWall/indirectWallPolyPatch.C
./meshes/polyMesh/polyPatches/derived/invisible/invisiblePolyPatch.C
./meshes/polyMesh/polyPatches/derived/inactive/inactivePolyPatch.C
./meshes/polyMesh/polyPatches/derived/inlet/inletPolyPatch.C
./meshes/polyMesh/polyPatches/derived/outlet/outletPolyPatch.C
./meshes/polyMesh/polyPatches/derived/nonConformalOrig/nonConformalOrigPolyPatch.C
./meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.C
./meshes/polyMesh/polyBoundaryMesh/polyBoundaryMeshEntries.C
./meshes/ProcessorTopology/commSchedule.C
./meshes/polyMesh/globalMeshData/globalMeshData.C
./meshes/polyMesh/globalMeshData/globalPoints.C
./meshes/polyMesh/globalMeshData/globalIndex.C
./meshes/polyMesh/syncTools/syncTools.C
./meshes/polyMesh/polyMeshTetDecomposition/polyMeshTetDecomposition.C
./meshes/polyMesh/polyMeshTetDecomposition/tetIndices.C
./meshes/polyMesh/zones/zone/zone.C
./meshes/polyMesh/zones/cellZone/cellZone.C
./meshes/polyMesh/zones/cellZone/cellZoneNew.C
./meshes/polyMesh/zones/faceZone/faceZone.C
./meshes/polyMesh/zones/faceZone/faceZoneNew.C
./meshes/polyMesh/zones/pointZone/pointZone.C
./meshes/polyMesh/zones/pointZone/pointZoneNew.C
./meshes/polyMesh/polyMesh.C
./meshes/polyMesh/polyMeshFromShapeMesh.C
./meshes/polyMesh/polyMeshIO.C
./meshes/polyMesh/polyMeshInitMesh.C
./meshes/polyMesh/polyMeshClear.C
./meshes/polyMesh/polyMeshUpdate.C
./meshes/polyMesh/polyMeshCheck/polyMeshCheck.C
./meshes/polyMesh/polyMeshCheck/polyMeshTools.C
./meshes/primitiveMesh/primitiveMesh.C
./meshes/primitiveMesh/primitiveMeshCellCells.C
./meshes/primitiveMesh/primitiveMeshCellCentresAndVols.C
./meshes/primitiveMesh/primitiveMeshCellEdges.C
./meshes/primitiveMesh/primitiveMeshCells.C
./meshes/primitiveMesh/primitiveMeshClear.C
./meshes/primitiveMesh/primitiveMeshEdgeCells.C
./meshes/primitiveMesh/primitiveMeshEdgeFaces.C
./meshes/primitiveMesh/primitiveMeshEdges.C
./meshes/primitiveMesh/primitiveMeshFaceCentresAndAreas.C
./meshes/primitiveMesh/primitiveMeshFindCell.C
./meshes/primitiveMesh/primitiveMeshPointCells.C
./meshes/primitiveMesh/primitiveMeshPointFaces.C
./meshes/primitiveMesh/primitiveMeshPointPoints.C
./meshes/primitiveMesh/primitiveMeshCellPoints.C
./meshes/primitiveMesh/primitiveMeshCalcCellShapes.C
./meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshCheck.C
./meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshCheckPointNearness.C
./meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshCheckEdgeLength.C
./meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.C
./meshes/primitiveMesh/primitivePatch/patchZones.C
./meshes/primitiveMesh/primitivePatch/walkPatch.C
./meshes/meshShapes/cellMatcher/cellMatcher.C
./meshes/meshShapes/cellMatcher/hexMatcher.C
./meshes/meshShapes/cellMatcher/wedgeMatcher.C
./meshes/meshShapes/cellMatcher/prismMatcher.C
./meshes/meshShapes/cellMatcher/pyrMatcher.C
./meshes/meshShapes/cellMatcher/tetWedgeMatcher.C
./meshes/meshShapes/cellMatcher/tetMatcher.C
./meshes/meshShapes/cellMatcher/degenerateMatcher.C
./meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.C
./meshes/polyMesh/polyTopoChangeMap/faceMapper/faceMapper.C
./meshes/polyMesh/polyTopoChangeMap/cellMapper/cellMapper.C
./meshes/polyMesh/polyDistributionMap/distributionMap.C
./meshes/polyMesh/polyDistributionMap/distributionMapBase.C
./meshes/polyMesh/polyDistributionMap/polyDistributionMap.C
./meshes/polyMesh/polyDistributionMap/IOdistributionMap.C
./meshes/polyMesh/polyTopoChangeMap/mapAddedPolyMesh.C
./meshes/polyMesh/polyMeshMap/polyMeshMap.C
./meshes/primitiveMesh/PrimitivePatch/PrimitivePatchName.C
./meshes/pointMesh/pointMesh.C
./meshes/pointMesh/pointMeshMapper/pointMapper.C
./meshes/pointMesh/pointMeshMapper/pointPatchMapper.C
./meshes/pointMesh/pointPatches/pointPatch/pointPatch.C
./meshes/pointMesh/pointPatches/facePointPatch/facePointPatch.C
./meshes/pointMesh/pointPatches/facePointPatch/facePointPatchNew.C
./meshes/pointMesh/pointPatches/basic/coupled/coupledPointPatch.C
./meshes/pointMesh/pointPatches/basic/generic/genericPointPatch.C
./meshes/pointMesh/pointPatches/constraint/cyclic/cyclicPointPatch.C
./meshes/pointMesh/pointPatches/constraint/cyclicSlip/cyclicSlipPointPatch.C
./meshes/pointMesh/pointPatches/constraint/empty/emptyPointPatch.C
./meshes/pointMesh/pointPatches/constraint/processor/processorPointPatch.C
./meshes/pointMesh/pointPatches/constraint/processorCyclic/processorCyclicPointPatch.C
./meshes/pointMesh/pointPatches/constraint/symmetryPlane/symmetryPlanePointPatch.C
./meshes/pointMesh/pointPatches/constraint/symmetry/symmetryPointPatch.C
./meshes/pointMesh/pointPatches/constraint/wedge/wedgePointPatch.C
./meshes/pointMesh/pointPatches/constraint/internal/internalPointPatch.C
./meshes/pointMesh/pointPatches/derived/coupled/coupledFacePointPatch.C
./meshes/pointMesh/pointPatches/derived/wall/wallPointPatch.C
./meshes/pointMesh/pointPatches/derived/indirectWall/indirectWallPointPatch.C
./meshes/pointMesh/pointPatches/derived/invisible/invisiblePointPatch.C
./meshes/pointMesh/pointPatches/derived/inactive/inactivePointPatch.C
./meshes/pointMesh/pointPatches/derived/inlet/inletPointPatch.C
./meshes/pointMesh/pointPatches/derived/outlet/outletPointPatch.C
./meshes/pointMesh/pointPatches/derived/nonConformalOrig/nonConformalOrigPointPatch.C
./meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.C
./meshes/boundBox/boundBox.C
./meshes/treeBoundBox/treeBoundBox.C
./meshes/meshTools/matchPoints.C
./meshes/meshTools/simpleVTKWriter.C
./fields/UniformDimensionedFields/uniformDimensionedFields.C
./fields/cloud/cloud.C
./fields/Fields/fieldMappers/directFieldMapper/directFieldMapper.C
./fields/Fields/fieldMappers/generalFieldMapper/generalFieldMapper.C
./fields/Fields/fieldMappers/setSizeFieldMapper/setSizeFieldMapper.C
./fields/Fields/labelField/labelField.C
./fields/Fields/scalarField/scalarField.C
./fields/Fields/sphericalTensorField/sphericalTensorField.C
./fields/Fields/diagTensorField/diagTensorField.C
./fields/Fields/symmTensorField/symmTensorField.C
./fields/Fields/tensorField/tensorField.C
./fields/Fields/quaternionField/quaternionField.C
./fields/Fields/triadField/triadField.C
./fields/Fields/complexFields/complexFields.C
./fields/Fields/labelField/labelIOField.C
./fields/Fields/labelField/labelFieldIOField.C
./fields/Fields/scalarField/scalarIOField.C
./fields/Fields/scalarField/scalarFieldIOField.C
./fields/Fields/vectorField/vectorIOField.C
./fields/Fields/vectorField/vectorFieldIOField.C
./fields/Fields/vector2DField/vector2DIOField.C
./fields/Fields/vector2DField/vector2DFieldIOField.C
./fields/Fields/sphericalTensorField/sphericalTensorIOField.C
./fields/Fields/sphericalTensorField/sphericalTensorFieldIOField.C
./fields/Fields/diagTensorField/diagTensorIOField.C
./fields/Fields/diagTensorField/diagTensorFieldIOField.C
./fields/Fields/symmTensorField/symmTensorIOField.C
./fields/Fields/symmTensorField/symmTensorFieldIOField.C
./fields/Fields/tensorField/tensorIOField.C
./fields/Fields/tensorField/tensorFieldIOField.C
./fields/Fields/quaternionField/quaternionIOField.C
./fields/Fields/triadField/triadIOField.C
./fields/Fields/transformField/transformField.C
./fields/Fields/symmTransformField/symmTransformFieldSpecialisations.C
./fields/pointPatchFields/basic/calculated/calculatedPointPatchFields.C
./fields/pointPatchFields/basic/coupled/coupledPointPatchFields.C
./fields/pointPatchFields/basic/value/valuePointPatchFields.C
./fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.C
./fields/pointPatchFields/basic/zeroGradient/zeroGradientPointPatchFields.C
./fields/pointPatchFields/constraint/cyclic/cyclicPointPatchFields.C
./fields/pointPatchFields/constraint/cyclicSlip/cyclicSlipPointPatchFields.C
./fields/pointPatchFields/constraint/empty/emptyPointPatchFields.C
./fields/pointPatchFields/constraint/processor/processorPointPatchFields.C
./fields/pointPatchFields/constraint/processorCyclic/processorCyclicPointPatchFields.C
./fields/pointPatchFields/constraint/symmetryPlane/symmetryPlanePointPatchFields.C
./fields/pointPatchFields/constraint/symmetry/symmetryPointPatchFields.C
./fields/pointPatchFields/constraint/wedge/wedgePointPatchFields.C
./fields/pointPatchFields/constraint/internal/internalPointPatchFields.C
./fields/pointPatchFields/derived/slip/slipPointPatchFields.C
./fields/pointPatchFields/derived/fixedNormalSlip/fixedNormalSlipPointPatchFields.C
./fields/pointPatchFields/derived/uniformFixedValue/uniformFixedValuePointPatchFields.C
./fields/pointPatchFields/derived/timeVaryingUniformFixedValue/timeVaryingUniformFixedValuePointPatchFields.C
./fields/pointPatchFields/derived/codedFixedValue/codedFixedValuePointPatchFields.C
./fields/GeometricFields/pointFields/pointFields.C
./fields/ReadFields/ReadFields.C
./meshes/bandCompression/bandCompression.C
./meshes/preservePatchTypes/preservePatchTypes.C
./interpolations/patchToPatchInterpolation/PatchToPatchInterpolationName.C
./interpolations/interpolationTable/tableReaders/tableReaders.C
./interpolations/interpolationTable/tableReaders/openFoam/openFoamTableReaders.C
./interpolations/interpolationTable/tableReaders/csv/csvTableReaders.C
./interpolations/interpolationWeights/interpolationWeights/interpolationWeights.C
./interpolations/interpolationWeights/linearInterpolationWeights/linearInterpolationWeights.C
./interpolations/interpolationWeights/splineInterpolationWeights/splineInterpolationWeights.C
./interpolations/primitivePatchInterpolation/primitivePatchInterpolation.C
./algorithms/indexedOctree/indexedOctreeName.C
./algorithms/indexedOctree/treeDataCell.C
./algorithms/indexedOctree/volumeType.C
./algorithms/dynamicIndexedOctree/dynamicIndexedOctreeName.C
./algorithms/dynamicIndexedOctree/dynamicTreeDataPoint.C
./algorithms/polygonTriangulate/polygonTriangulate.C
./algorithms/serialThreads/serialThreads.C
./primitives/BlockCoeff/blockCoeffBase.C
./primitives/BlockCoeff/scalarBlockCoeff.C
./primitives/BlockCoeff/sphericalTensorBlockCoeff.C
./primitives/BlockCoeff/symmTensorBlockCoeff.C
./primitives/BlockCoeff/tensorBlockCoeff.C
./fields/expandContract/expandTensorField.C
./fields/CoeffField/scalarCoeffField.C
./fields/CoeffField/sphericalTensorCoeffField.C
./fields/CoeffField/symmTensorCoeffField.C
./fields/CoeffField/tensorCoeffField.C
./matrices/blockLduMatrix/BlockLduMatrix/blockLduMatrices.C
./matrices/blockLduMatrix/BlockLduMatrix/scalarBlockLduMatrix.C
./matrices/blockLduMatrix/BlockLduMatrix/sphericalTensorBlockLduMatrix.C
./matrices/blockLduMatrix/BlockLduMatrix/symmTensorBlockLduMatrix.C
./matrices/blockLduMatrix/BlockLduMatrix/tensorBlockLduMatrix.C
./matrices/blockLduMatrix/BlockLduMatrix/BlockConstraint/scalarBlockConstraint.C
./matrices/blockLduMatrix/BlockLduPrecons/BlockLduPrecon/blockLduPrecons.C
./matrices/blockLduMatrix/BlockLduPrecons/BlockNoPrecon/blockNoPrecons.C
./matrices/blockLduMatrix/BlockLduSmoothers/BlockLduSmoother/blockLduSmoothers.C
./matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/scalarBlockSolverPerformance.C
./matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/blockLduSolvers.C
./matrices/blockLduMatrix/BlockLduSolvers/BlockDiagonal/blockDiagonalSolvers.C
./matrices/blockLduMatrix/BlockLduMatrix/BlockLduInterface/BlockLduInterfaceFields.C
./coordinate/systems/frameData.C
./coordinate/systems/coordinateSystem.C
./coordinate/systems/coordinateSystemNew.C
./coordinate/systems/coordinateSystemTransform.C
./coordinate/systems/coordinateSystems.C
./coordinate/systems/cartesianCS.C
./coordinate/systems/cylindricalCS.C
./coordinate/systems/sphericalCS.C
./coordinate/rotation/axesRotation.C
./coordinate/rotation/axisAngleRotation.C
./coordinate/rotation/coordinateRotation.C
./coordinate/rotation/cylindricalRotation.C
./coordinate/rotation/identityRotation.C
./coordinate/rotation/EulerCoordinateRotation.C
./coordinate/rotation/STARCDCoordinateRotation.C
./VectorN/primitives/vector1/diagTensor1.C
./VectorN/primitives/vector1/sphericalTensor1.C
./VectorN/primitives/vector1/tensor1.C
./VectorN/primitives/vector1/vector1.C
./VectorN/primitives/vector2/diagTensor2.C
./VectorN/primitives/vector2/sphericalTensor2.C
./VectorN/primitives/vector2/tensor2.C
./VectorN/primitives/vector2/vector2.C
./VectorN/primitives/vector4/diagTensor4.C
./VectorN/primitives/vector4/sphericalTensor4.C
./VectorN/primitives/vector4/tensor4.C
./VectorN/primitives/vector4/vector4.C
./VectorN/primitives/vector5/diagTensor5.C
./VectorN/primitives/vector5/sphericalTensor5.C
./VectorN/primitives/vector5/tensor5.C
./VectorN/primitives/vector5/vector5.C
./VectorN/primitives/vector6/diagTensor6.C
./VectorN/primitives/vector6/sphericalTensor6.C
./VectorN/primitives/vector6/tensor6.C
./VectorN/primitives/vector6/vector6.C
./VectorN/primitives/vector8/diagTensor8.C
./VectorN/primitives/vector8/sphericalTensor8.C
./VectorN/primitives/vector8/tensor8.C
./VectorN/primitives/vector8/vector8.C
./VectorN/Lists/VectorNLists.C
./VectorN/DimensionedTypes/dimensionedVectorTensorN.C
./VectorN/Fields/VectorNFields.C
./VectorN/Fields/TensorNFields.C
./VectorN/Fields/DiagTensorNFields.C
./VectorN/Fields/SphericalTensorNFields.C
./VectorN/expandContract/ExpandTensorNField.C
./VectorN/BlockLduMatrices/blockVectorNSolvers.C
./VectorN/BlockLduMatrices/blockVectorNMatrices.C
./VectorN/BlockLduMatrices/blockVectorNLduInterfaceFields.C
./graph/curve/curve.C
./graph/graph.C
./graph/writers/rawGraph/rawGraph.C
./graph/writers/gnuplotGraph/gnuplotGraph.C
./graph/writers/xmgrGraph/xmgrGraph.C
./graph/writers/jplotGraph/jplotGraph.C
./meshes/data/data.C
$<TARGET_OBJECTS:OSspecific>
)

# Append to list of Helyx targets for emake autocompletion
get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
list(APPEND TARGETS ${TARGET_NAME})
set_property(GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS ${TARGETS})

add_helyx_include_target(${CMAKE_CURRENT_LIST_DIR})
add_to_dir_target(${TARGET_NAME})

target_include_directories(${TARGET_NAME}
    PUBLIC $<BUILD_INTERFACE:${HELYX_PROJECT_DIR}/src/OSspecific/${HELYX_SYSTEM_NAME}>
    PUBLIC $<INSTALL_INTERFACE:${HELYX_DIR_NAME}/src/OSspecific/${HELYX_SYSTEM_NAME}>
    PUBLIC ${THIRDPARTY_MPI_INC}
)


helyx_link_libraries(${TARGET_NAME}
    PRIVATE ${THIRDPARTY_ZLIB}
    PUBLIC ${THIRDPARTY_MPI}
)

if("${HELYX_SYSTEM_NAME}" STREQUAL "POSIX")
    # OpenFOAM and OSspecific directories are added by default, so
    # call the function here, but don't need to specify any directories
    helyx_additional_includes(${TARGET_NAME}
        PRIVATE
            ${THIRDPARTY_ZLIB_INC}
    )
    target_sources(${TARGET_NAME}
        PRIVATE ./primitives/random/random11/random11.C
        PRIVATE ./global/profiling/profiling.C
        PRIVATE ./global/profiling/profilingInformation.C
        PRIVATE ./global/profiling/profilingSysInfo.C
        PRIVATE ./global/profiling/profilingTrigger.C
    )

    # Put libOpenFOAM.so in the mpi subdirectory, to allow users to store
    # compilations against several different MPI's at once
    set_target_properties(${TARGET_NAME}
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}
    )

    target_link_options(${TARGET_NAME}
        PRIVATE
            -lrt -ldl -lm
    )

    if (BFD_LIBRARY)
        target_link_libraries(${TARGET_NAME} PRIVATE ${BFD_LIBRARY})
    endif()

elseif("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")

    helyx_link_libraries(${TARGET_NAME}
        PUBLIC ${THIRDPARTY_REGEX}
        PUBLIC ${THIRDPARTY_STACK_TRACE}
    )

    # This is horrible, but it's what Allwmake does and it seems to work.
    # At some stage, I'll get round to doing something nicer...
    add_custom_command(TARGET ${TARGET_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/flexInclude/"
        COMMAND ${CMAKE_COMMAND} -E copy "/usr/include/FlexLexer.h" "${CMAKE_CURRENT_BINARY_DIR}/flexInclude/"
        COMMENT "Copying FlexLexer.h to OpenFOAM/flexInclude"
    )

    # Put libOpenFOAM.so in the mpi subdirectory, to allow users to store
    # compilations against several different MPI's at once
    # CMake treats Windows .dll files as runtime targets
    set_target_properties(${TARGET_NAME}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}
    )

    target_include_directories("${TARGET_NAME}"
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/flexInclude/>
    )
endif()


get_directory_property(OF_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_CURRENT_LIST_DIR} COMPILE_DEFINITIONS)
set(OF_COMPILE_DEFINITIONS "-D${OF_COMPILE_DEFINITIONS}")
string (REPLACE ";" " -D" OF_COMPILE_DEFINITIONS "${OF_COMPILE_DEFINITIONS}")

string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)


# Configure the dynamicCode Makefile template.
# It requires OF_COMPILE_DEFINITIONS, so, this must be done here!
# Cannot use $MPI_ARCH_PATH/include because we could in fact have $MPI_ARCH_PATH/**/include.
# So, make the include paths relative to $MPI_ARCH_PATH.
# Also, need to parse THIRDPARTY_MPI_INC if multiple paths are found.
set(DYNCODE_THIRDPARTY_MPI_INC "${THIRDPARTY_MPI_INC}")
string (REPLACE "${MPI_ARCH_PATH}" "\$(MPI_ARCH_PATH)"
    DYNCODE_THIRDPARTY_MPI_INC "${DYNCODE_THIRDPARTY_MPI_INC}"
)
string (REPLACE ";" " -I"
    DYNCODE_THIRDPARTY_MPI_INC "${DYNCODE_THIRDPARTY_MPI_INC}"
)
configure_file(
    ${HELYX_PROJECT_DIR}/etc/cmake/templates/template-dynamicCode.Makefile_
    ${HELYX_PROJECT_DIR}/etc/codeTemplates/dynamicCode/makefiles/${HELYX_OPTIONS}/dynamicCode.Makefile
    NEWLINE_STYLE UNIX
)

if("${TBB_FOUND}")
    helyx_link_libraries(${TARGET_NAME}
        PUBLIC ${THIRDPARTY_TBB}
    )
endif()

#if("${OPENVDB_FOUND}")
#    target_sources(${TARGET_NAME}
#        PRIVATE ./algorithms/VDBTree/VDBTree.C
#    )
#
#    helyx_additional_includes(${TARGET_NAME}
#        ${THIRDPARTY_OPENVDB_INC}
#    )
#
#    target_link_libraries(${TARGET_NAME}
#        PUBLIC ${THIRDPARTY_OPENVDB}
#    )
#
#    set_target_properties(
#        ${TARGET_NAME}
#        PROPERTIES
#        CXX_STANDARD 14
#    )
#endif()


# Manually override installation instructions here - need to change output directory
# For Windows, only pack runtime (.exe and .dll, NOT .dll.a)
if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    set(RUNTIME_TOGGLE RUNTIME)
endif()

# This is a core component
install(TARGETS ${TARGET_NAME}
    EXPORT HELYX-core
    ${RUNTIME_TOGGLE}
        DESTINATION ${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}/platforms/${HELYX_OPTIONS}/lib/${HELYX_MPI_NAME}
        COMPONENT ${TARGET_PLATFORM_NAME}
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
)

# Now that OpenFOAM lives in the $HELYX_MPI_NAME subdir, we need to add
# $ORIGIN/.. and $ORIGIN to the RUNPATH in order for run-time selection to work
get_target_property(install_runpath ${TARGET_NAME} INSTALL_RPATH)
string(REPLACE "$ORIGIN/../lib"
    "$ORIGIN;$ORIGIN/.." install_runpath
    "${install_runpath}"
    )
set_target_properties(${TARGET_NAME}
    PROPERTIES INSTALL_RPATH "${install_runpath}"
    )


# Append to list of Helyx targets for emake autocompletion
get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
list(APPEND TARGETS ${TARGET_NAME})
set_property(GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS ${TARGETS})
