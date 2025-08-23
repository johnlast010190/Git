#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2019 Engys Ltd

[----------------------------------------------------------------------------]]

set(TARGET_NAME "swakFunctionObjects")

add_helyx_library(${TARGET_NAME} SHARED
./manipulateField.C
./manipulateFieldFunctionObject.C
./manipulateTopoSet.C
./manipulateTopoSetFunctionObject.C
./manipulatePatchField.C
./manipulatePatchFieldFunctionObject.C
./expressionField.C
./expressionFieldFunctionObject.C
./clearExpressionField.C
./clearExpressionFieldFunctionObject.C
./createSampledSurface.C
./createSampledSurfaceFunctionObject.C
./createSampledSet.C
./createSampledSetFunctionObject.C
./addGlobalVariable.C
./addGlobalVariableFunctionObject.C
./removeGlobalVariable.C
./removeGlobalVariableFunctionObject.C
./swakCodedFunctionObject.C
./calculateGlobalVariables.C
./calculateGlobalVariablesFunctionObject.C
./solverPerformanceToGlobalVariables.C
./solverPerformanceToGlobalVariablesFunctionObject.C
./addForeignMeshes.C
./addForeignMeshesFunctionObject.C
./foreignMeshesFollowTime.C
./foreignMeshesFollowTimeFunctionObject.C
./solvePDECommon.C
./solvePDECommonFiniteVolume.C
./solveLaplacianPDE.C
./solveLaplacianPDEFunctionObject.C
./solveTransportPDE.C
./solveTransportPDEFunctionObject.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        surfMesh-includes
        fileFormats-includes
        ${HELYX_PROJECT_DIR}/modules/swak4Foam/Libraries/filesThatAreOnlyInSomeDistributions/OutputFilerDropin
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE fvOptions
    PRIVATE sampling
    PRIVATE utilityFunctionObjects
    PRIVATE dynamicFvMesh
    PRIVATE topoChangerFvMesh
    PRIVATE dynamicMesh
    PRIVATE lagrangian
)
helyx_link_libraries_with_generated_headers(${TARGET_NAME} PRIVATE swak4FoamParsers)
