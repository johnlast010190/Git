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
    (c) 2019-2022 Engys Ltd

[----------------------------------------------------------------------------]]

set(TARGET_NAME "swak4FoamParsers")

helyx_flex_target(FieldValueExpressionLexer
    ./FieldValueExpressionLexer.ll
    ${CMAKE_CURRENT_BINARY_DIR}/FieldValueExpressionLexer.C
    COMPILE_FLAGS "-f"
)

helyx_flex_target(PatchValueExpressionLexer
    ./PatchValueExpressionLexer.ll
    ${CMAKE_CURRENT_BINARY_DIR}/PatchValueExpressionLexer.C
    COMPILE_FLAGS "-f"
)

helyx_flex_target(SubsetValueExpressionLexer
    ./SubsetValueExpressionLexer.ll
    ${CMAKE_CURRENT_BINARY_DIR}/SubsetValueExpressionLexer.C
    COMPILE_FLAGS "-f"
)

add_helyx_library(${TARGET_NAME} SHARED
    ../filesThatAreOnlyInSomeDistributions/OutputFilerDropin/outputFilterOutputControl.C
    ./helpers/printSwakVersion.C
    ./helpers/exprString.C
    ./helpers/switchableNotImplemented.C
    ./helpers/DebugOStream.C
    ./entryToExpression/entryToExpression.C
    ./entryToExpression/directEntryToExpression.C
    ./entryToExpression/stringEntryToExpression.C
    ./entryToExpression/boolEntryToExpression.C
    ./entryToExpression/vectorEntryToExpression.C
    ./entryToExpression/tensorEntryToExpression.C
    ./entryToExpression/symmTensorEntryToExpression.C
    ./entryToExpression/sphericalTensorEntryToExpression.C
    ./entryToExpression/dimensionedScalarEntryToExpression.C
    ./entryToExpression/dimensionedVectorEntryToExpression.C
    ./entryToExpression/dimensionedTensorEntryToExpression.C
    ./entryToExpression/dimensionedSymmTensorEntryToExpression.C
    ./entryToExpression/dimensionedSphericalTensorEntryToExpression.C
    ./repositories/RepositoryBase.C
    ./repositories/SurfacesRepository.C
    ./repositories/MeshesRepository.C
    ./repositories/SetsRepository.C
    ./repositories/GlobalVariablesRepository.C
    ./namedEnums/MeshInterpolationOrder.C
    ./namedEnums/LogicalAccumulationNamedEnum.C
    ./namedEnums/NumericAccumulationNamedEnum.C
    ./ExpressionDriverWriter.C
    ./ExpressionResult/ExpressionResult.C
    ./ExpressionResult/StoredExpressionResult.C
    ./ExpressionResult/DelayedExpressionResult.C
    ./ExpressionResult/StackExpressionResult.C
    ./ExpressionResult/StoredStackExpressionResult.C
    ./setsAndSurfaces/swakRegistryProxySurface.C
    ./setsAndSurfaces/swakRegistryProxySet.C
    ./makeSwakDataEntry.C
    ./CommonValueExpressionDriver.C
    ${FLEX_FieldValueExpressionLexer_OUTPUTS}
    ./FieldValueExpressionDriver.C
    ${FLEX_PatchValueExpressionLexer_OUTPUTS}
    ./PatchValueExpressionDriver.C
    ${FLEX_SubsetValueExpressionLexer_OUTPUTS}
    ./SubsetValueExpressionDriver.C
    ./SetSubsetValueExpressionDriver.C
    ./CellSetValueExpressionDriver.C
    ./FaceSetValueExpressionDriver.C
    ./CellZoneValueExpressionDriver.C
    ./FaceZoneValueExpressionDriver.C
    ./SampledSurfaceValueExpressionDriver.C
    ./SampledSetValueExpressionDriver.C
    ./initSwakFunctionObject.C
    ./plugins/CommonPluginFunction.C
    ./plugins/FieldValuePluginFunction.C
    ./plugins/PatchValuePluginFunction.C
    ./plugins/CellSetValuePluginFunction.C
    ./plugins/CellZoneValuePluginFunction.C
    ./plugins/FaceSetValuePluginFunction.C
    ./plugins/FaceZoneValuePluginFunction.C
    ./plugins/SampledSetValuePluginFunction.C
    ./plugins/SampledSurfaceValuePluginFunction.C
)

target_add_bison_sources(${TARGET_NAME}
    FieldValueExpressionParser
    PatchValueExpressionParser
    SubsetValueExpressionParser
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fileFormats-includes
        ${HELYX_PROJECT_DIR}/modules/swak4Foam/Libraries/filesThatAreOnlyInSomeDistributions/OutputFilerDropin
    PUBLIC
        lagrangian-includes
        lagrangianIntermediate-includes
        fluidThermophysicalModels-includes
)

# Let CMake deals with Flex headers unless FLEX_ARCH_PATH is set.
# Make sure the right Flex headers are included in case FLEX_ARCH_PATH is set
# and a different Flex version exist in /usr/*
#
# N.B.: Looks like there is a bug in CMake find Flex module when using tool-chains/windows:
# Flex executable is found and the right includes were used but
# FLEX_INCLUDE_DIR is reported as not found (FLEX_INCLUDE_DIR-NOTFOUND)
# Therefore, only include the header if FLEX_ARCH_PATH is set for linux
if(NOT "" STREQUAL "${FLEX_ARCH_PATH}" AND NOT "${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    helyx_additional_includes(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_FLEX_INC}
    )
endif()

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PUBLIC sampling
    PRIVATE lagrangian
    PRIVATE dynamicMesh
    PRIVATE surfMesh
    PRIVATE meshTools
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/FaceZoneValueExpressionDriver.C
    ${CMAKE_CURRENT_LIST_DIR}/CellZoneValueExpressionDriver.C
    ${CMAKE_CURRENT_LIST_DIR}/namedEnums/LogicalAccumulationNamedEnum.C
    ${CMAKE_CURRENT_LIST_DIR}/namedEnums/NumericAccumulationNamedEnum.C
    ${CMAKE_CURRENT_LIST_DIR}/FaceSetValueExpressionDriver.C
    ${BISON_FieldValueExpressionParser_OUTPUT_SOURCE}
    ${BISON_PatchValueExpressionParser_OUTPUT_SOURCE}
    ${BISON_SubsetValueExpressionParser_OUTPUT_SOURCE}
    ${CMAKE_CURRENT_LIST_DIR}/CommonValueExpressionDriver.C
    ${FLEX_FieldValueExpressionLexer_OUTPUTS}
    ${CMAKE_CURRENT_LIST_DIR}/FieldValueExpressionDriver.C
    ${FLEX_PatchValueExpressionLexer_OUTPUTS}
    ${CMAKE_CURRENT_LIST_DIR}/PatchValueExpressionDriver.C
    ${FLEX_SubsetValueExpressionLexer_OUTPUTS}
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
