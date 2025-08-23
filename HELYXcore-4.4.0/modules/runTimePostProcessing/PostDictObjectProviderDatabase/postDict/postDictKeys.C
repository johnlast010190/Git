/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
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
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

namespace Foam::functionObjects::runTimeVis
{

namespace genericKeys
{
    const char* NAME_KEY = "name";
    const char* VISIBLE_KEY = "visible";
    const char* MIN_KEY = "min";
    const char* MAX_KEY = "max";
    const char* POSITION_KEY = "position";
    const char* RADIUS_KEY = "radius";
    const char* CENTER_KEY = "center";
    const char* ORIGIN_KEY = "origin";
    const char* AXIS_KEY = "axis";
    const char* NORMAL_KEY = "normal";
    const char* TYPE_KEY = "type";
    const char* POINT_1_KEY = "point1";
    const char* POINT_2_KEY = "point2";
    const char* WIDGETS_COLOUR_KEY = "widgetsColor";

    const char* UNIFORM_FOLDER_NAME = "uniform";
    const char* SCALAR_SUBDICT_KEY = "scalar";
}

namespace idKeys
{

    const char* NAME_KEY = genericKeys::NAME_KEY;
    const char* REGION_KEY = "region";
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
}

namespace postDictKeys
{
    const char* VISUALISATIONS_SUBDICT_KEY = "visualisations";
    const char* OPERATING_POINT_SUBDICT_KEY = "operatingPointDict";
    const char* INSTANCE_KEY = "instance";
    const char* OBJECTS_SUBDICT_KEY = "objects";
    const char* OBJECTS_TYPE_KEY = genericKeys::TYPE_KEY;
    const char* GROUPS_SUBDICT_KEY = "groups";
    const char* REFERENCE_FRAMES_SUBDICT_KEY = "referenceFrames";
}

namespace operatingPointKeys
{
    const char* PARENT_MESH_NAME_KEY = "parentMeshName";
}

namespace controlDictKeys
{
    const char* FUNCTIONS_SUBDICT_KEY = "functions";
    const char* FUNCTIONS_TYPE_KEY = genericKeys::TYPE_KEY;
    const char* EXPORT_FORMATS_KEY = "exportFormats";
    const char* DEBUG_KEY = "debug";
}

namespace geometryDictKeys
{
    const char* SURFACES_SUBDICT_KEY = "surfaces";
    const char* LINES_SUBDICT_KEY = "lines";
    const char* FORMAT_KEY = "format";
}

namespace imageExportKeys
{
    const char* TRANSPARENT_KEY = "transparentBackground";
    const char* CROP_TO_CONTENTS_KEY = "cropToContents";
    const char* CENTER_CROPPED_CONTENTS_KEY = "centerCroppedContents";
    const char* IMAGE_WIDTH_KEY = "width";
    const char* IMAGE_HEIGHT_KEY = "height";
    const char* EDF_FIELDS_KEY = "edfFields";
    const char* EDF_COMPRESSION_LEVEL_KEY = "edfCompressionLevel";
}

namespace pvdExportKeys
{
    const char* PVD_FIELDS_TYPE_KEY = "pvdFieldsType";
    const char* PVD_FIELDS_LIST_KEY = "pvdFieldsList";
    const char* PVD_WRITE_POINT_DATA_KEY = "writePointData";
    const char* PVD_COMPRESSION_LEVEL_KEY = "pvdCompressionLevel";
    const char* PVD_WRITE_ONE_FILE_PER_PROCESS_KEY = "pvdOneFilePerProcess";
    const char* PVD_REMOVE_OLD_OBSOLETE_FILES_KEY = "pvdRemoveObsoleteFiles";
}

namespace visualisationKeys
{
    const char* GEOMETRIES_SUBDICT_KEY = "geometry";
    const char* REGIONS_SUBDICT_KEY = "regions";
    const char* OBJECTS_SUBDICT_KEY = postDictKeys::OBJECTS_SUBDICT_KEY;
    const char* GROUPS_SUBDICT_KEY = postDictKeys::GROUPS_SUBDICT_KEY;
    const char* REPRESENTATION_KEY = "representation";
    const char* VISIBLE_KEY = "visible";
    const char* OPACITY_KEY = "opacity";
    const char* COLOR_FIELD_KEY = "colorField";
    const char* COLOR_KEY = "color";
    const char* LINE_COLOR_KEY = "lineColor";
    const char* LINE_THICKNESS_KEY = "lineThickness";
    const char* BACKFACE_STYLING_KEY = "backfaceStyling";
    const char* SHOW_ACTIVE_KEY = "showActiveColor";
    const char* ACTIVE_COLOR_KEY = "activeColor";
    const char* SHOW_NORMALS_KEY = "showNormals";
    const char* NORMALS_LENGTH_KEY = "normalsLength";
    const char* NORMALS_RATIO_KEY = "normalsRatio";
    const char* NORMALS_COLOR_KEY = "normalsColor";
    const char* NORMALS_OPACITY_KEY = "normalsOpacity";
}
namespace surfaceLICKeys
{
    const char* SHOW_KEY = "showSurfaceLIC";
    const char* STEP_SIZE_KEY = "stepSize";
    const char* NUMBER_OF_STEPS_KEY = "numberOfSteps";
    const char* VECTOR_FIELD_KEY = "vectorArray";
}


namespace regionKeys
{
    const char* FACE_ZONES_KEY = "faceZones";
    const char* FILE_SOURCES_KEY = "fileSources";
    const char* CELL_ZONES_KEY = "cellZones";
    const char* INTERNAL_BOUNDARIES_KEY = "internalBoundaries";
    const char* DEFAULT_REGION_KEY = "region0";
    const char* VOLUME_MESH_NAME = "volumeMesh";
}

namespace boundaryGroupKeys
{
    const char* EXTERNAL_BOUNDARY_GROUP_NAME = "externalBoundaries";
    const char* INTERNAL_BOUNDARY_GROUP_NAME = "internalBoundaries";
}

namespace colourMapKeys
{
    const char* COLOR_SPACE_KEY = "colorSpace";
    const char* HSV_WRAP_KEY = "hsvWrap";
    const char* POINTS_KEY = "points";
    const char* COLORS_KEY = "colors";
}

namespace colourLegendsKeys
{
    const char* TITLE_KEY = "title";
    const char* THICKNESS_KEY = "thickness";
    const char* HEIGHT_KEY = "height";
    const char* LOCATION_KEY = genericKeys::POSITION_KEY;
    const char* COORDINATES_KEY = "coordinates";
    const char* VERTICAL_KEY = "vertical";
    const char* LABEL_FORMAT_KEY = "labelFormat";
    const char* NUMBER_OF_LABELS_KEY = "numberOfLabels";
    const char* AUTOMATIC_LABELS_KEY = "automaticLabels";
    const char* SHOW_TICKS_KEY = "showTicks";
    const char* FONT_DICT_KEY = "legendLabelFont";
}

namespace fontKeys {

    const char* FONT_KEY = "font";
    const char* SIZE_KEY = "fontSize";
    const char* COLOR_KEY = "textColor";
    const char* OPACITY_KEY = "textOpacity";
    const char* BOLD_KEY = "bold";
    const char* ITALICS_KEY = "italics";
    const char* SHADOW_KEY = "shadow";
}

namespace colourLookupTableKeys
{
    const char* COLOUR_LOOKUP_DICT_KEY = "colorMap";
    const char* TYPE_KEY = "colorLegendType";
    const char* RESOLUTION_KEY = "resolution";
    const char* INVERTED_KEY = "inverted";
    const char* AUTOMATIC_RANGE_KEY = "automaticRange";
    const char* RANGE_KEY = "range";
}

namespace surfaceKeys
{
    const char* MIN_KEY = genericKeys::MIN_KEY;
    const char* MAX_KEY = genericKeys::MAX_KEY;
    const char* POINT_1_KEY = genericKeys::POINT_1_KEY;
    const char* POINT_2_KEY = genericKeys::POINT_2_KEY;
    const char* RADIUS_KEY = genericKeys::RADIUS_KEY;
    const char* PARENT_NAME_KEY = "parentName";
    const char* POINT_AND_NORMAL_DICT_KEY = "pointAndNormalDict";
    const char* BASE_POINT_KEY = "basePoint";
    const char* NORMAL_VECTOR_KEY = "normalVector";
    const char* DIAGONAL_KEY = "diagonal";
    const char* INNER_RADIUS_KEY = "innerRadius";
    const char* OUTER_RADIUS_KEY = "outerRadius";
    const char* CENTRE_KEY = "centre";
    const char* HEX_MESH_GEOMETRY_SUBDICT_KEY = visualisationKeys::GEOMETRIES_SUBDICT_KEY;
    const char* CASTELLATED_MESH_CONTROLS_SUBDICT_KEY = "castellatedMeshControls";
    const char* FEATURES_LIST_KEY = "features";
    const char* FILE_KEY = "file";
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
    const char* SOLIDS_KEY = "solids";
    const char* NAME_KEY = idKeys::NAME_KEY;
    const char* SURFACE_TYPE_SUBDICT_KEY = "surfaceType";
    const char* MERGE_COPLANAR_KEY = "mergeCoplanar";
}

namespace surfaceTransformationKeys
{
    const char* TRANSFORMS_KEY = "transforms";
    const char* TYPE_KEY = "type";
    const char* ORIGIN_KEY = "aboutPoint";
    const char* SCALE_KEY = "scaleVec";
    const char* TRANSLATION_KEY = "translateVec";
    const char* ROTATIONS_KEY = "rollPitchYaw";
}

namespace blockMeshDictKeys
{
    const char* VERTICES_KEY = "vertices";
    const char* BLOCKS_KEY = "blocks";
    const char* PATCHES_KEY = "patches";
}

namespace cuttingTypeKeys
{
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
    const char* POSITION_KEY = genericKeys::POSITION_KEY;
    const char* ROTATION_KEY = "rotation";
    const char* SCALE_KEY = "scale";
    const char* CENTER_KEY = genericKeys::CENTER_KEY;
    const char* RADIUS_KEY = genericKeys::RADIUS_KEY;
    const char* AXIS_KEY = genericKeys::AXIS_KEY;
    const char* STL_FILE_KEY = "stlFile";
    const char* ORIGIN_KEY = genericKeys::ORIGIN_KEY;
    const char* NORMAL_KEY = genericKeys::NORMAL_KEY;
}

namespace sliceObjectKeys
{
    const char* CRINKLE_KEY = "crinkle";
    const char* SLICE_TYPE_DICT_KEY = "sliceType";
    const char* OFFSETS_KEY = "offsets";
}

namespace clipObjectKeys
{
    const char* CRINKLE_KEY = sliceObjectKeys::CRINKLE_KEY;
    const char* INSIDE_OUT_KEY = "insideOut";
    const char* SLICE_TYPE_DICT_KEY = sliceObjectKeys::SLICE_TYPE_DICT_KEY;
}

namespace streamlinesSourceKeys
{
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
    const char* SOURCE_DICT_KEY = "seed";
    const char* NUMBER_OF_POINTS_KEY = "numberOfPoints";
    const char* POINT_1_KEY = genericKeys::POINT_1_KEY;
    const char* POINT_2_KEY = genericKeys::POINT_2_KEY;
    const char* PATCH_ID_KEY = "patchId";
    const char* RADIUS_KEY = genericKeys::RADIUS_KEY;
    const char* CENTER_KEY = genericKeys::CENTER_KEY;
}

namespace streamlinesObjectKeys
{
    const char* VECTOR_FIELD_KEY = "vectorField";
    const char* MAX_LENGTH_KEY = "maxLength";
    const char* MAX_STEPS_KEY = "maxSteps";
    const char* RADIUS_KEY = genericKeys::RADIUS_KEY;
}

namespace fieldSamplingObjectKeys
{
    const char* SAMPLING_BOUNDS_MIN_KEY = "samplingMin";
    const char* SAMPLING_BOUNDS_MAX_KEY = "samplingMax";
    const char* ELEMENTS_KEY = "elements";
    const char* FIELDS_KEY = "fields";
    const char* TOLERANCE_TYPE = "toleranceType";
    const char* USE_TOLERANCE = "useTolerance";
    const char* TOLERANCE_VALUE = "toleranceValue";
}

namespace line3DObjectKeys
{
    const char* SOURCE_TYPE_DATA_KEY = "line3dType";
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
    const char* N_SAMPLING_POINTS_KEY = "nSamplingPoints";
    const char* POINT_1_KEY = genericKeys::POINT_1_KEY;
    const char* POINT_2_KEY = genericKeys::POINT_2_KEY;
    const char* SURFACE_INTERSECTION_KEY = "surfaceIntersection";
    const char* FEATURE_LINE_ID_KEY = "featureLineId";
    const char* FILE_NAME_KEY = genericKeys::NAME_KEY;
}

namespace sceneKeys
{
    const char* SCENES_DICT_KEY = "scenes";
    const char* BACKGROUND_COLORS_DICT_KEY = "colors";
    const char* WIDGETS_DICT_KEY = "colors";
    const char* CAMERAS_DICT_KEY = "cameras";
    const char* ITEMS_DICT_KEY = "items";
    const char* COLOR_LEGEND_OPTIONS_DICT_KEY = "colorLegendOptions";
    const char* COLOR_LEGENDS_DICT_KEY = "colorLegends";
}

namespace fieldSamplingFunctionObjectKeys
{
    const char* SAMPLING_FILE_KEY = "samplingFile";
}

namespace backgroundKeys
{
    const char* BACKGROUND_1_KEY = "background1";
    const char* BACKGROUND_2_KEY = "background2";
}

namespace cameraKeys
{
    const char* FOCAL_POINT_KEY = "focalPoint";
    const char* POSITION_KEY = genericKeys::POSITION_KEY;
    const char* PARALLEL_PROJECTION_KEY = "parallelProjection";
    const char* UP_KEY = "up";
    const char* PARALLEL_SCALE_KEY = "parallelScale";
    const char* NAME_KEY = idKeys::NAME_KEY;
    const char* COORDINATE_TYPE_KEY = "coordinateType";
    const char* REFERENCE_FRAME_KEY = "referenceFrame";
}

namespace axisWidgetKeys
{
    const char* COLOUR_KEY = genericKeys::WIDGETS_COLOUR_KEY;
    const char* VISIBLE_KEY = "showAxes";
}

namespace gridWidgetKeys
{
    const char* GRID_WIDGET_OPTIONS_DICT_KEY = "cubeAxesOptions";
    const char* SHOW_AXES_GRID_KEY         = "showAxesGrid";
    const char* AXES_GRID_COLOR_KEY        = "axesGridColor";
    const char* AXES_GRID_LABEL_COLOR_KEY  = "axesGridLabelColor";
    const char* SHOW_AXES_GRID_TITLE_KEY   = "showAxesGridTitle";
    const char* X_AXIS_TITLE_KEY           = "XAxisTitle";
    const char* Y_AXIS_TITLE_KEY           = "YAxisTitle";
    const char* Z_AXIS_TITLE_KEY           = "ZAxisTitle";
    const char* MANUAL_FACE_SELECT_KEY     = "manualFaceSelect";
    const char* SHOW_XYMIN_KEY             = "showXYmin";
    const char* SHOW_XYMAX_KEY             = "showXYmax";
    const char* SHOW_YZMIN_KEY             = "showYZmin";
    const char* SHOW_YZMAX_KEY             = "showYZmax";
    const char* SHOW_ZXMIN_KEY             = "showZXmin";
    const char* SHOW_ZXMAX_KEY             = "showZXmax";
    const char* USE_CUSTOM_BOUNDS_KEY      = "useCustomBounds";
    const char* X_BOUNDS_KEY               = "XBounds";
    const char* Y_BOUNDS_KEY               = "YBounds";
    const char* Z_BOUNDS_KEY               = "ZBounds";
    const char* USE_CUSTOM_LABELS_KEY      = "useCustomLabels";
    const char* USE_X_LABELS_KEY           = "useXDistanceLabel";
    const char* USE_Y_LABELS_KEY           = "useYDistanceLabel";
    const char* USE_Z_LABELS_KEY           = "useZDistanceLabel";
    const char* X_DISTANCE_LABEL_KEY       = "useXDistance";
    const char* Y_DISTANCE_LABEL_KEY       = "useYDistance";
    const char* Z_DISTANCE_LABEL_KEY       = "useZDistance";
}

namespace logoWidgetKeys
{
    const char* VISIBLE_KEY = "showLogo";
}

namespace timestepWidgetKeys
{
    const char* COLOUR_KEY = genericKeys::WIDGETS_COLOUR_KEY;
    const char* VISIBLE_KEY = "showTimestep";
    const char* FORMAT_KEY = "timestepLabelFormat";
    const char* COORDINATES_KEY = "timestepLabelCoordinates";
    const char* FONT_DICT_KEY = "timestepLabelFont";
}

namespace vectorWidgetKeys
{
    const char* VECTOR_WIDGET_OPTIONS_DICT_KEY = "vectorWidget";

    const char* VISIBLE_KEY = "showVectorWidget";
    const char* VECTOR_REPRESENTATION_KEY = "vector";
    const char* VECTOR_3D_REP_OPTIONS_KEY = "vector3DRepOptions";
    const char* VECTOR_TEXT_OPTIONS_KEY = "vectorTextRepOptions";

    const char* SHOW_VECTOR_3D_REP_KEY = "showVector3DRep";
    const char* SHOW_VECTOR_RESULTANT_KEY = "showResultant";
    const char* SHOW_VECTOR_COMPONENTS_KEY = "showComponents";
    const char* VECTOR_WINDOW_SIZE_KEY = "scale";

    const char* SHOW_TITLE_KEY = "showTitle";
    const char* COORDINATES_KEY = "coordinates";

    const char* SHOW_VECTOR_TEXT_KEY = "showVectorTextRep";
    const char* LABEL_FORMAT_KEY = "labelFormat";
    const char* FONT_SIZE_KEY = "fontSize";
    const char* SHOW_TEXT_BOX_KEY = "showTextBox";

    const char* RUNTIME_INFO_DICT = "runTimeInfo";
    const char* RESULTS_SUBDICT = "results";
    const char* INTERFOAM_SUBDICT = "interFoam";
    const char* FLOWSOLVER_SUBDICT = "flowSolver";
    const char* MAX_FRAME_ACCELERATION_COMPONENTS_KEY = "maxFrameAccelerationComponents";
    const char* MIN_FRAME_ACCELERATION_COMPONENTS_KEY = "minFrameAccelerationComponents";
    const char* FRAME_ACCELERATION_KEY = "frameAcceleration";
    const char* GRAVITY_ACCELERATION_KEY = "g";

    const char* VECTOR_SUBDICT_KEY = "vector";
}

namespace referenceFramesWidgetKeys
{
const char *REFERENCE_FRAMES_WIDGET_LIST_KEY = "referenceFrames";
}

namespace activePatchKeys
{
    const char* ACTIVE_FIELD_KEY = "activeFieldData";
}

namespace meshObjectsKeys
{
    const char* TYPE_KEY = genericKeys::TYPE_KEY;
    const char* REFERENCE_FRAME_ORIGIN_KEY = genericKeys::ORIGIN_KEY;
    const char* REFERENCE_FRAME_E1_KEY = "e1";
    const char* REFERENCE_FRAME_E2_KEY = "e2";
    const char* REFERENCE_FRAME_E3_KEY = "e3";
    const char* REFERENCE_FRAME_TYPE_KEY = "coordinateFrame";
    const char* MOTION_COORDINATE_FRAME_TYPE_KEY = "motionCoordinateFrame";
    const char* COORDINATE_SYSTEM_KEY = "coordinateSystem";
    const char* GLOBAL_COORDINATE_SYSTEM_KEY = "Global";
    const char* INVALID_COORDINATE_SYSTEM_KEY = "Invalid";
}

namespace turbopropKeys
{
    const char* INLET_PATCHES_KEY = "inletPatches";
    const char* OUTLET_PATCHES_KEY = "outletPatches";
    const char* HUB_PATCHES_KEY = "hubPatches";
    const char* SHROUD_PATCHES_KEY = "shroudPatches";
    const char* STREAM_ELEMENTS_KEY = "streamElements";
    const char* SPAN_ELEMENTS_KEY = "spanElements";
    const char* ORIGIN_KEY = genericKeys::ORIGIN_KEY;
    const char* AXIS_KEY = genericKeys::AXIS_KEY;
    const char* X_DIRECTION_KEY = "xDirection";
    const char* SPAN_VALUE_KEY = "spanValue";
    const char* STREAM_VALUE_KEY = "streamValue";
    const char* BLADE_TO_BLADE_KEY = "bladeToBlade";
    const char* BLADE_PATCHES_KEY = "bladeLoadingPatches";
}

namespace groupKeys
{
    const char* SOURCES_KEY = "sources";
}

namespace thresholdKeys
{
    const char* FIELD_KEY = "field";
    const char* ALL_POINTS_CRITERION_KEY = "allPointsCriterion";
    const char* MIN_THRESHOLD_KEY = "minThreshold";
    const char* MAX_THRESHOLD_KEY = "maxThreshold";
    const char* INVERT_KEY = "invert";
}

namespace transformKeys
{
    const char* TRANSFORMATIONS_LIST_KEY = "transformations";
    const char* TRANSFORMATION_TYPE_KEY = "transformationType";
    const char* INCLUDE_SOURCE_KEY = "includeSource";
    const char* COPY_COUNT_KEY = "copyCount";
    const char* OFFSET_DEGREES_KEY = "offsetDegrees";
    const char* ORIGIN_KEY = genericKeys::ORIGIN_KEY;
    const char* AXIS_KEY = genericKeys::AXIS_KEY;
    const char* OFFSET_KEY = "offset";
    const char* SLICE_TYPE_DICT_KEY = sliceObjectKeys::SLICE_TYPE_DICT_KEY;
    const char* NORMAL_KEY = genericKeys::NORMAL_KEY;
    const char* DISTANCE_KEY = "distance";
    const char* ANGLE_DEGREES_KEY = "angleDegrees";
    const char* DIRECTION_KEY = "direction";
    const char* RATIO_KEY = "ratio";
}

namespace importDataObjectKeys
{
    const char* PATH_KEY = "path";
    const char* FORMAT_KEY = "format";
    const char* TIME_STEP_KEY = "timestep";
    const char* PVD_KEY = "pvd";
    const char* VTK_KEY = "vtk";
    const char* FIELD_SAMPLING_KEY = "fieldSampling";
    const char* ENSIGHT_KEY = "ensight";

}

namespace objectKeys
{
    const char* SOURCE_KEY = "source";
}

namespace tool3Dto2D
{
    const char* SOURCES_KEY = groupKeys::SOURCES_KEY;
    const char* ORIGIN_KEY = genericKeys::ORIGIN_KEY;
    const char* NORMAL_KEY = genericKeys::NORMAL_KEY;
}

namespace surfaceSplitterKeys
{
    const char* ANGLE_KEY = transformKeys::ANGLE_DEGREES_KEY;
}


} // End namespace Foam
