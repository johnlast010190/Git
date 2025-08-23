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
    (c) 2019-2025 Engys Ltd

[----------------------------------------------------------------------------]]

set(TARGET_NAME "simpleFunctionObjects")

if (${HELYX_SYSTEM_NAME} STREQUAL "POSIX")
    set(POSIX_ONLY_SRC
    ./misc/writeOldTimesOnSignal/writeOldTimesOnSignal.C
    ./misc/provokeSignal/provokeSignalFunctionObject.C
    )
else()
    message(STATUS "Not compiling SWAK signal function objects for Windows")
endif()

add_helyx_library(${TARGET_NAME} SHARED
./SolverInfo/SolverInfo.C
./DataStructures/TimeClone/TimeClone.C
./DataStructures/TimeClone/TimeCloneList.C
./DataStructures/TimelineCollection/TimelineCollection.C
./simpleFunctionObject/simpleFunctionObject.C
./simpleDataFunctionObject/simpleDataFunctionObject.C
./timelineFunctionObject/timelineFunctionObject.C
./distributionFunctionObject/distributionFunctionObject.C
./distributionFunctionObject/field/fieldDistributionFunctionObject.C
./distributionFunctionObject/patchField/patchFieldDistributionFunctionObject.C
./patch/patchFunctionObject/patchFunctionObject.C
./patch/patchFieldFunctionObject/patchFieldFunctionObject.C
./patch/patchAverageFunctionObject/patchAverageFunctionObject.C
./patch/patchIntegrateFunctionObject/patchIntegrateFunctionObject.C
./patch/patchMassFlowFunctionObject/patchMassFlowFunctionObject.C
./patch/patchMassFlowAverageFunctionObject/patchMassFlowAverageFunctionObject.C
./patch/patchFieldFlowFunctionObject/patchFieldFlowFunctionObject.C
./volume/volumeAverageFunctionObject/volumeAverageFunctionObject.C
./volume/volumeFieldFunctionObject/volumeFieldFunctionObject.C
./volume/volumeIntegrateFunctionObject/volumeIntegrateFunctionObject.C
./volume/volumeMinMaxFunctionObject/volumeMinMaxFunctionObject.C
./misc/panicDump/panicDumpFunctionObject.C
./misc/listRegisteredObjects/listRegisteredObjectsFunctionObject.C
./misc/listMeshData/listMeshDataFunctionObject.C
./misc/readGravitation/readGravitationFunctionObject.C
./misc/readAndUpdateFields/readAndUpdateFields.C
./misc/readAndUpdateFields/readAndUpdateFieldsFunctionObject.C
${POSIX_ONLY_SRC}
./misc/loadTopoSets/loadTopoSetsFunctionObject.C
./writeFields/writeAdditionalFields/writeAdditionalFieldsFunctionObject.C
./writeFields/writeFieldsOften/writeFieldsOftenFunctionObject.C
./writeFields/writeFieldsGeneral/writeFieldsGeneralFunctionObject.C
./writeFields/conditionDrivenWriting/conditionDrivenWritingFunctionObject.C
./writeFields/conditionDrivenWriting/writeIfFieldOutside/writeIfFieldOutsideFunctionObject.C
./trackDictionary/trackDictionaryFunctionObject.C
./functionObjectProxy/functionObjectListProxy/functionObjectListProxy.C
./functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.C
./functionObjectProxy/executeIfExecutableFitsFunctionObject/executeIfExecutableFitsFunctionObject.C
./functionObjectProxy/executeIfObjectExistsFunctionObject/executeIfObjectExistsFunctionObject.C
./functionObjectProxy/executeIfEnvironmentVariableFunctionObject/executeIfEnvironmentVariableFunctionObject.C
./functionObjectProxy/executeIfFunctionObjectPresentFunctionObject/executeIfFunctionObjectPresentFunctionObject.C
./functionObjectProxy/executeIfParallelSerialFunctionObject/executeIfParallelSerialFunctionObject.C
./functionObjectProxy/executeIfStartTimeFunctionObject/executeIfStartTimeFunctionObject.C
./functionObjectProxy/executeIfWriteTimeFunctionObject/executeIfWriteTimeFunctionObject.C
./functionObjectProxy/executeIfOpenFOAMVersionBiggerEqualFunctionObject/executeIfOpenFOAMVersionBiggerEqualFunctionObject.C
./functionObjectProxy/executeIfSetOrZoneFunctionObject/executeIfSetOrZoneFunctionObject.C
./functionObjectProxy/executeIfPatchFunctionObject/executeIfPatchFunctionObject.C
./functionObjectProxy/dynamicFunctionObjectListProxy/dynamicFunctionObjectListProxy.C
./functionObjectProxy/dynamicFunctionObjectListProxy/providers/fromFileDictionaryProvider.C
./functionObjectProxy/dynamicFunctionObjectListProxy/providers/stdoutFromCommandProvider.C
./misc/writeAndEndFunctionObject/writeAndEndFunctionObject.C
./misc/writeAndEndFieldRangeFunctionObject/writeAndEndFieldRangeFunctionObject.C
./modelLoading/loadThermoModelFunctionObjects.C
./modelLoading/loadIncompressibleTurbulenceModelFunctionObject.C
./modelLoading/loadCompressibleTurbulenceModelFunctionObject.C
./modelLoading/loadSLGThermoModelFunctionObject.C
./timeManipulation/timeManipulationFunctionObject/timeManipulationFunctionObject.C
./timeManipulation/setDeltaTByTimelineFunctionObject/setDeltaTByTimelineFunctionObject.C
./manipulateFvSolutionFvSchemes/manipulateFvSolutionFvSchemesFunctionObject/manipulateFvSolutionFvSchemesFunctionObject.C
./manipulateFvSolutionFvSchemes/manipulateFvSolutionFvSchemesFunctionObject/manipulateFvSolutionFvSchemesFunctionObjectRereadMethods.C
./manipulateFvSolutionFvSchemes/timeDependentFvSolutionFvSchemesFunctionObject/timeDependentFvSolutionFvSchemesFunctionObject.C
./initialization/initPotentialFlow/initPotentialFlowFunctionObject.C
./update/updateSimple/updateSimpleFunctionObject.C
./update/recalcPhi/recalcPhiFunctionObject.C
./update/correctThermo/correctThermoFunctionObject.C
./update/recalcThermoH/recalcThermoHFunctionObject.C
./update/recalcThermoHs/recalcThermoHsFunctionObject.C
./update/recalcThermoE/recalcThermoEFunctionObject.C
./update/recalcThermoHe/recalcThermoHeFunctionObject.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        surfMesh-includes
        fluidThermophysicalModels-includes
        thermophysicalProperties-includes
        radiationModels-includes
        regionModels-includes
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
        ${HELYX_PROJECT_DIR}/src/transportModels
        ${HELYX_PROJECT_DIR}/modules/swak4Foam/Libraries/filesThatAreOnlyInSomeDistributions/OutputFilerDropin
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE sampling
    PRIVATE fluidThermophysicalModels
    PRIVATE incompressibleTransportModels
    PRIVATE incompressibleTurbulenceModels
    PRIVATE turbulenceModels
    PRIVATE compressibleTurbulenceModels
    PRIVATE turbulenceModels
    PRIVATE specie
    PRIVATE lagrangian
    PRIVATE lagrangianIntermediate
    PRIVATE SLGThermo
)
helyx_link_libraries_with_generated_headers(${TARGET_NAME} PRIVATE swak4FoamParsers)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/manipulateFvSolutionFvSchemes/manipulateFvSolutionFvSchemesFunctionObject/manipulateFvSolutionFvSchemesFunctionObjectRereadMethods.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
