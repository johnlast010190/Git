#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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

# If it reaches this file, is because runTimePostProcessing_FOUND is TRUE

# RTPP fails to compile without EVTK, so EVTK_REQUIRED should default to ON
# unless otherwise set (note that this cache variable doesn't use FORCE).
set(VTK_REQUIRED ON CACHE STRING "")

# EVTK is mandatory for runTimePostProcessing compilation
# so, first set some sensible variables based on possible user combinations of
# ${module}_REQUIRED and ${dependency}_REQUIRED
helyx_check_module_mandatory_dependency(runTimePostProcessing VTK)

# if runTimePostProcessing compilation is disabled, simply return from this script
if(runTimePostProcessing_REQUIRED STREQUAL "OFF" OR NOT runTimePostProcessing_FOUND)
    return()
endif()

# check whether EVTK is available and set appropriate messages
helyx_add_module_dependency(runTimePostProcessing VTK)

# Return before processing all the CMakeLists
if(NOT VTK_FOUND AND NOT VTK_REQUIRED STREQUAL "OFF")
    return()
endif()

# ------------------------------------------------------------------------------

message(STATUS "runTimePostProcessing:")
message(CLEAN "    Using EVTK ext/ from ${VTK_DIR}")

# Append EVTK libs to RPATH, for both OpenGL and OSMesa backends
append_library_paths_to_rpaths(VTK)
strip_system_paths_from_thirdparty_includes(VTK)

# Also append LLVM and Mesa to RPATH for OSMesa backend
if("OSMesa" STREQUAL "${VTK_RENDERING_BACKEND}")
    #file(RELATIVE_PATH OSMesa_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${VTK_DIR}/EVTK-OSMesa)
    #file(RELATIVE_PATH OSMesaX11_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${VTK_DIR}/EVTK-OSMesa/X11)
    #list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/${OSMesaX11_dir}" "$ORIGIN/${OSMesa_dir}")
    file(RELATIVE_PATH LLVM_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${VTK_DIR}/LLVM/lib)
    file(RELATIVE_PATH MESA_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${VTK_DIR}/Mesa/lib)
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/${LLVM_dir}" "$ORIGIN/${MESA_dir}")
endif()

# Need to remove duplicates from CMAKE_INSTALL_RPATH
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

# ------------------------------------------------------------------------------

set(RTPP_POST_DICT_OBJECTS_DIR "${CMAKE_CURRENT_LIST_DIR}/PostDictObjectProviderDatabase")
set(RTPP_VISUALISATION_DIR "${CMAKE_CURRENT_LIST_DIR}/RunTimeVisualisation")
set(SURFACE_STATISTICS_DIR "${CMAKE_CURRENT_LIST_DIR}/SurfaceStatistics")
set(FIELD_SAMPLING_DIR "${CMAKE_CURRENT_LIST_DIR}/FieldSampling")
set(ANIMATIONS_DIR "${CMAKE_CURRENT_LIST_DIR}/Animations")

# First, generate availableColourmaps.H when necessary, with correct
# dependencies
add_custom_command(
    OUTPUT ${RTPP_VISUALISATION_DIR}/colourMaps/availableColourmaps.H
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/CMakeGenAvailableColourmaps.cmake
    DEPENDS
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blackBlueWhite
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blackBodyRadiation
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blackOrangeWhite
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blackToWhite
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToRedDES
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToRedDIV
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToRedHSV
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToRedJet
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToRedRGB
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToYellowDIV
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToYellowHSV
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/blueToYellowRGB
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/coldAndHot
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/coolToWarm
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/CpX
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/engys
        ${RTPP_VISUALISATION_DIR}/colourMaps/definitions/viridis
    COMMENT "Generating availableColourmaps.H"
)


set(TARGET_NAME "runTimeVisualisation")

add_helyx_library(${TARGET_NAME} SHARED
    ${CMAKE_CURRENT_LIST_DIR}/rtppFunctionObject.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/foamField.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/foamFields.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/meshAndFields.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/itemRequirements.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/visualisation/visualisation.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/baseClasses/visualisation/surfaceLIC.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/widgets/colourLegendData.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/geometry/blockMesh/blockMeshStream.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/scene/cameraData.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/referenceFrame/referenceFrame.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/geometry/surfaceTransformations/surfaceScaleData.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/geometry/surfaceTransformations/surfaceTranslateData.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/geometry/surfaceTransformations/surfaceRotateData.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/dataStructs/geometry/surfaceTransformationsData.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/dictionaries/dictionaries.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/fileSource/pvdFileDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/boundarySurfaceDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/featureLineDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceBoxDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceCylinderDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceFromFileDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfacePlaneDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceRingDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceSphereDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/geometry/surfaceDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/group/groupDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/cellZoneDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/faceZoneDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/meshDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/patchDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/processorBoundaryDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/mesh/volMeshDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/missingItem/missingItemDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/cutting/cuttingPlaneProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/cutting/cuttingBoxProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/cutting/cuttingCylinderProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/cutting/cuttingFileProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/cutting/cuttingSphereProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/line3d/line3DDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/line3d/featureLine3DDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/line3d/fromFileLine3DDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/line3d/straightLine3DDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/line3d/surfaceIntersectionLine3DDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/streamlineSources/lineSourceProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/streamlineSources/patchSourceProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/streamlineSources/pointCloudSourceProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/clipDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/fieldSamplingDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/glyphDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/isosurfaceDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/sliceDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/streamlinesDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/thresholdDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/transformDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/meridionalGridDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/turboSliceStreamwiseDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/turboSliceSpanwiseDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/bladeLoadingLineDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/inletToOutletLineDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/turboPost/hubToShroudLineDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/objects/importDataDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/testObjects/surfaceSplitterDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/externalItemDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/itemDataSetProvider.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/itemDataSetProviders/testObjects/tool3Dto2DDataSetProvider.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/postDict/postDictKeys.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/itemInfo.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/itemInfoFactory.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/baseItemInfoReader.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/objects/cutting/cuttingInfoFactory.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/objects/line3DSources/line3DSourceInfoFactory.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/objects/line3DSources/surfaceIntersectionLine3DSourceInfo.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/items/objects/streamlineSources/streamlineSourceInfoFactory.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/colourLookupTablesInfo.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/sceneInfo.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/infos/sceneInfos.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/storage/foamMeshes.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/storage/externalFields.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/storage/itemStorage.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/storage/itemInfoIndex.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/storage/referenceFrames.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/types/surfaceFormat.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/uniformDistribution/boundsBasedUniformInput.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/uniformDistribution/gatheredDomainUniformInput.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/uniformDistribution/spreadDomainUniformInput.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/Utils.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/ParallelUtils.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/patchTypes.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/vtkDataHandlingTools.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/boundsUtils.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/basicProfiler.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/Utils/memoryMonitor.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/vtkWrappers/engysLabelCoreDataArray.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/vtkWrappers/engysScalarCoreDataArray.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/vtkWrappers/engysVectorCoreDataArray.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/vtkWrappers/rtppVtkVtuAdaptor.C

    ${RTPP_POST_DICT_OBJECTS_DIR}/databaseInstance/postDictObjectProviderDatabaseInstance.C
      ${RTPP_POST_DICT_OBJECTS_DIR}/databaseInstance/instanceServices.C
    ${RTPP_POST_DICT_OBJECTS_DIR}/postDictObjectProviderDatabase.C


    # RunTimeVisualisation
    ${RTPP_VISUALISATION_DIR}/colourLookupTable/colourLookupTable.C
    ${RTPP_VISUALISATION_DIR}/colourLookupTable/colourLookupTableProvider.C

    ${RTPP_VISUALISATION_DIR}/colourMaps/colours.C

    ${RTPP_VISUALISATION_DIR}/exporters/pvdWriter.C

    ${RTPP_VISUALISATION_DIR}/infos/pvdInfo.C
    ${RTPP_VISUALISATION_DIR}/infos/renderInfo.C
    ${RTPP_VISUALISATION_DIR}/infos/runTimeVisualisationInfo.C

    ${RTPP_VISUALISATION_DIR}/rendering/compositers/allGatherTransparentCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/basicCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/edfCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/iceTCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/masterGatherCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/compositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/compositerFactory.C
    ${RTPP_VISUALISATION_DIR}/rendering/compositers/parallelDepthPeelCompositer.C
    ${RTPP_VISUALISATION_DIR}/rendering/edfRenderers/edfWriter.C
    ${RTPP_VISUALISATION_DIR}/rendering/pngRenderers/pngWriter.C
    ${RTPP_VISUALISATION_DIR}/rendering/renderManager.C
    ${RTPP_VISUALISATION_DIR}/rendering/rtppActor.C

    ${RTPP_VISUALISATION_DIR}/scene/scene.C
    ${RTPP_VISUALISATION_DIR}/scene/item.C
    ${RTPP_VISUALISATION_DIR}/scene/items.C

    ${RTPP_VISUALISATION_DIR}/widgets/widgetsInOverlay.C
    ${RTPP_VISUALISATION_DIR}/widgets/widgetsInScene.C
    ${RTPP_VISUALISATION_DIR}/widgets/inScene/gridWidget.C
    ${RTPP_VISUALISATION_DIR}/widgets/inScene/referenceFrameWidget.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/axesWidget.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/colorLegend.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/colorLegends.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/logoWidget.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/timestepWidget.C
    ${RTPP_VISUALISATION_DIR}/widgets/overlay/vectorWidget.C

    ${RTPP_VISUALISATION_DIR}/colourMaps/availableColourmaps.H

    ${RTPP_VISUALISATION_DIR}/runTimeVisualisation.C


    # Surface Statistics
    ${SURFACE_STATISTICS_DIR}/surfaceStatistics.C

    ${SURFACE_STATISTICS_DIR}/surfaceStatisticsWriter.C
    ${SURFACE_STATISTICS_DIR}/controlDict/controlDictKeys.C

    # Field Sampling
    ${FIELD_SAMPLING_DIR}/fieldSampling.C

    ${FIELD_SAMPLING_DIR}/fieldSamplingWriter.C
    ${FIELD_SAMPLING_DIR}/controlDict/controlDictKeysFS.C

    # Animations
    ${ANIMATIONS_DIR}/animation.C
    ${ANIMATIONS_DIR}/controlDict/controlDictKeysAN.C
)

# We can only build RTPP against he Engys VTK. So I guess this is not required anymore
#if(NOT ${USING_ENGYS_VTK})
#    message(STATUS "Linking runTimeVisualisation functionObject against non-standard VTK")
#    include(${VTK_USE_FILE})
#endif()

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        OSspecific-includes
)

target_include_directories(${TARGET_NAME}
    SYSTEM PRIVATE ${THIRDPARTY_VTK_INC}
)

target_include_directories(${TARGET_NAME}
    PRIVATE ${RTPP_POST_DICT_OBJECTS_DIR}
    PRIVATE ${RTPP_VISUALISATION_DIR}
)

target_link_libraries(${TARGET_NAME}
    PRIVATE surfMesh
    PRIVATE dynamicMesh
    PRIVATE conversion
    PRIVATE finiteVolume
    ${THIRDPARTY_VTK}
)

# MemoryMonitor requires psapi in Windows
string( TOLOWER "${HELYX_SYSTEM_NAME}" HELYX_SYSTEM_NAME_LOWER )
if ("mswindows" STREQUAL "${HELYX_SYSTEM_NAME_LOWER}")
    target_link_libraries(${TARGET_NAME}
            PRIVATE psapi.lib
    )
endif()

# RTPP can no longer build against these older versions. So I guess this is not required anymore
## Specialise for VTK version 7 to improve forwards-compatability
#if(${VTK_VERSION_8})
#    target_compile_definitions(${TARGET_NAME} PRIVATE VTK_VERSION_8)
#endif()


# Use RPATH instead of RUNPATH, because RPATH is inherited by inherited
# dependencies, whereas RUNPATH is not.  This allows the VTK libraries that
# this target links to to use this target's RPATH variables to find OSMesa
# libraries.
if("POSIX" STREQUAL ${HELYX_SYSTEM_NAME})
    set_target_properties(${TARGET_NAME}
        PROPERTIES LINK_FLAGS -Wl,--disable-new-dtags
    )
endif()


# Add install rpath to point to standard GUI library location
get_property(RTPP_INSTALL_RPATH
    TARGET ${TARGET_NAME}
    PROPERTY INSTALL_RPATH
)

if (DEFINED ENV{HELYX_EXT_PATH})
    set(EXT_PATH "$ENV{HELYX_EXT_PATH}")
    message(CLEAN "    Also adding the system EVTK ext/ path to rpath")
else()
    set(EXT_PATH "$ORIGIN/../../../../../GUI/ext")
    message(CLEAN "    Also adding the hardcoded EVTK ext/ path to rpath")
endif()

set(RTPP_INSTALL_RPATH
    "${EXT_PATH}/EVTK-${VTK_RENDERING_BACKEND};${RTPP_INSTALL_RPATH}"
)
if("OSMesa" STREQUAL "${VTK_RENDERING_BACKEND}")
    # This shouldn't be necessary, because the OSMesa archpath should be set
    # properly when EVTK is compiled
    set(RTPP_INSTALL_RPATH
        "${EXT_PATH}/EVTK-${VTK_RENDERING_BACKEND}/X11;${RTPP_INSTALL_RPATH}"
    )
    set(RTPP_INSTALL_RPATH
        "${EXT_PATH}/LLVM/lib;${RTPP_INSTALL_RPATH}"
    )
    set(RTPP_INSTALL_RPATH
        "${EXT_PATH}/Mesa/lib;${RTPP_INSTALL_RPATH}"
    )
endif()

set_target_properties(${TARGET_NAME}
    PROPERTIES INSTALL_RPATH "${RTPP_INSTALL_RPATH}"
)
