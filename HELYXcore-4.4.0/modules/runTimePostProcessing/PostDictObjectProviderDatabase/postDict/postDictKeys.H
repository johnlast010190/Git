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
    (c) 2019 OpenCFD Ltd.
    (c) 2023-2025 Engys Ltd.

Namespace
    Foam::functionObjects::runTimeVis::postDictKeys

Description
    Contains the definition of the keys used in the postDict

SourceFiles
    <postDictKeys.C>

\*---------------------------------------------------------------------------*/

#pragma once

namespace Foam::functionObjects::runTimeVis
{

namespace idKeys
{
    extern const char* NAME_KEY;
    extern const char* REGION_KEY;
    extern const char* TYPE_KEY;
}

namespace postDictKeys
{
    extern const char* VISUALISATIONS_SUBDICT_KEY;
    extern const char* OPERATING_POINT_SUBDICT_KEY;
    extern const char* INSTANCE_KEY;
    extern const char* OBJECTS_SUBDICT_KEY;
    extern const char* OBJECTS_TYPE_KEY;
    extern const char* GROUPS_SUBDICT_KEY;
    extern const char* REFERENCE_FRAMES_SUBDICT_KEY;
}

namespace operatingPointKeys
{
    extern const char* PARENT_MESH_NAME_KEY;
}

namespace controlDictKeys
{
    extern const char* FUNCTIONS_SUBDICT_KEY;
    extern const char* FUNCTIONS_TYPE_KEY;
    extern const char* EXPORT_FORMATS_KEY;
    extern const char* DEBUG_KEY;
}

namespace geometryDictKeys
{
    extern const char* SURFACES_SUBDICT_KEY;
    extern const char* LINES_SUBDICT_KEY;
    extern const char* FORMAT_KEY;
}

namespace imageExportKeys
{
    extern const char* TRANSPARENT_KEY;
    extern const char* CROP_TO_CONTENTS_KEY;
    extern const char* CENTER_CROPPED_CONTENTS_KEY;
    extern const char* IMAGE_WIDTH_KEY;
    extern const char* IMAGE_HEIGHT_KEY;
    extern const char* EDF_FIELDS_KEY;
    extern const char* EDF_COMPRESSION_LEVEL_KEY;
}

namespace pvdExportKeys
{
    extern const char* PVD_FIELDS_TYPE_KEY;
    extern const char* PVD_FIELDS_LIST_KEY;
    extern const char* PVD_WRITE_POINT_DATA_KEY;
    extern const char* PVD_COMPRESSION_LEVEL_KEY;
    extern const char* PVD_WRITE_ONE_FILE_PER_PROCESS_KEY;
    extern const char* PVD_REMOVE_OLD_OBSOLETE_FILES_KEY;
}

namespace visualisationKeys
{
    extern const char* GEOMETRIES_SUBDICT_KEY;
    extern const char* REGIONS_SUBDICT_KEY;
    extern const char* OBJECTS_SUBDICT_KEY;
    extern const char* GROUPS_SUBDICT_KEY;
    extern const char* REPRESENTATION_KEY;
    extern const char* VISIBLE_KEY;
    extern const char* OPACITY_KEY;
    extern const char* COLOR_FIELD_KEY;
    extern const char* COLOR_KEY;
    extern const char* LINE_COLOR_KEY;
    extern const char* LINE_THICKNESS_KEY;
    extern const char* BACKFACE_STYLING_KEY;
    extern const char* SHOW_ACTIVE_KEY;
    extern const char* ACTIVE_COLOR_KEY;
    extern const char* SHOW_NORMALS_KEY;
    extern const char* NORMALS_LENGTH_KEY;
    extern const char* NORMALS_RATIO_KEY;
    extern const char* NORMALS_COLOR_KEY;
    extern const char* NORMALS_OPACITY_KEY;
}

namespace surfaceLICKeys
{
    extern const char* SHOW_KEY;
    extern const char* STEP_SIZE_KEY;
    extern const char* NUMBER_OF_STEPS_KEY;
    extern const char* VECTOR_FIELD_KEY;
}

namespace regionKeys
{
    extern const char* FACE_ZONES_KEY;
    extern const char* FILE_SOURCES_KEY;
    extern const char* CELL_ZONES_KEY;
    extern const char* INTERNAL_BOUNDARIES_KEY;
    extern const char* DEFAULT_REGION_KEY;
    extern const char* VOLUME_MESH_NAME;
}

namespace boundaryGroupKeys
{
    extern const char* EXTERNAL_BOUNDARY_GROUP_NAME;
    extern const char* INTERNAL_BOUNDARY_GROUP_NAME;
}

namespace colourMapKeys
{
    extern const char* COLOR_SPACE_KEY;
    extern const char* HSV_WRAP_KEY;
    extern const char* POINTS_KEY;
    extern const char* COLORS_KEY;
}

namespace colourLegendsKeys
{
    extern const char* TITLE_KEY;
    extern const char* THICKNESS_KEY;
    extern const char* HEIGHT_KEY;
    extern const char* LOCATION_KEY;
    extern const char* COORDINATES_KEY;
    extern const char* VERTICAL_KEY;
    extern const char* LABEL_FORMAT_KEY;
    extern const char* NUMBER_OF_LABELS_KEY;
    extern const char* AUTOMATIC_LABELS_KEY;
    extern const char* SHOW_TICKS_KEY;
    extern const char* FONT_DICT_KEY;
}

namespace fontKeys {

    extern const char* FONT_KEY;
    extern const char* SIZE_KEY;
    extern const char* COLOR_KEY;
    extern const char* OPACITY_KEY;
    extern const char* BOLD_KEY;
    extern const char* ITALICS_KEY;
    extern const char* SHADOW_KEY;
}

namespace colourLookupTableKeys
{
    extern const char* COLOUR_LOOKUP_DICT_KEY; // Inside the colour legend dicts
    extern const char* TYPE_KEY;
    extern const char* RESOLUTION_KEY;
    extern const char* INVERTED_KEY;
    extern const char* AUTOMATIC_RANGE_KEY;
    extern const char* RANGE_KEY;
}

namespace surfaceKeys
{
    extern const char* MIN_KEY;
    extern const char* MAX_KEY;
    extern const char* POINT_1_KEY;
    extern const char* POINT_2_KEY;
    extern const char* RADIUS_KEY;
    extern const char* PARENT_NAME_KEY;
    extern const char* POINT_AND_NORMAL_DICT_KEY;
    extern const char* BASE_POINT_KEY;
    extern const char* NORMAL_VECTOR_KEY;
    extern const char* DIAGONAL_KEY;
    extern const char* INNER_RADIUS_KEY;
    extern const char* OUTER_RADIUS_KEY;
    extern const char* CENTRE_KEY;
    extern const char* HEX_MESH_GEOMETRY_SUBDICT_KEY;
    extern const char* CASTELLATED_MESH_CONTROLS_SUBDICT_KEY;
    extern const char* FEATURES_LIST_KEY;
    extern const char* FILE_KEY;
    extern const char* TYPE_KEY;
    extern const char* SOLIDS_KEY;
    extern const char* NAME_KEY;
    extern const char* SURFACE_TYPE_SUBDICT_KEY;
    extern const char* MERGE_COPLANAR_KEY;
}

namespace surfaceTransformationKeys
{
    extern const char* TRANSFORMS_KEY;
    extern const char* TYPE_KEY;
    extern const char* ORIGIN_KEY;
    extern const char* SCALE_KEY;
    extern const char* TRANSLATION_KEY;
    extern const char* ROTATIONS_KEY;
}

namespace blockMeshDictKeys
{
    extern const char* VERTICES_KEY;
    extern const char* BLOCKS_KEY;
    extern const char* PATCHES_KEY;
}

namespace cuttingTypeKeys
{
    extern const char* TYPE_KEY;
    extern const char* POSITION_KEY;
    extern const char* ROTATION_KEY;
    extern const char* SCALE_KEY;
    extern const char* CENTER_KEY;
    extern const char* RADIUS_KEY;
    extern const char* AXIS_KEY;
    extern const char* STL_FILE_KEY;
    extern const char* ORIGIN_KEY;
    extern const char* NORMAL_KEY;
}

namespace sliceObjectKeys
{
    extern const char* CRINKLE_KEY;
    extern const char* OFFSETS_KEY;
    extern const char* SLICE_TYPE_DICT_KEY;
}

namespace clipObjectKeys
{
    extern const char* CRINKLE_KEY;
    extern const char* INSIDE_OUT_KEY;
    extern const char* SLICE_TYPE_DICT_KEY;
}

namespace streamlinesSourceKeys
{
    extern const char* TYPE_KEY;
    extern const char* SOURCE_DICT_KEY;
    extern const char* NUMBER_OF_POINTS_KEY;
    extern const char* POINT_1_KEY;
    extern const char* POINT_2_KEY;
    extern const char* PATCH_ID_KEY;
    extern const char* RADIUS_KEY;
    extern const char* CENTER_KEY;
}

namespace streamlinesObjectKeys
{
    extern const char* VECTOR_FIELD_KEY;
    extern const char* MAX_LENGTH_KEY;
    extern const char* MAX_STEPS_KEY;
    extern const char* RADIUS_KEY;
}

namespace fieldSamplingObjectKeys
{
    extern const char* SAMPLING_BOUNDS_MIN_KEY;
    extern const char* SAMPLING_BOUNDS_MAX_KEY;
    extern const char* ELEMENTS_KEY;
    extern const char* FIELDS_KEY;
    extern const char* TOLERANCE_TYPE;
    extern const char* USE_TOLERANCE;
    extern const char* TOLERANCE_VALUE;
}

namespace line3DObjectKeys
{
    extern const char* SOURCE_TYPE_DATA_KEY;
    extern const char* TYPE_KEY;
    extern const char* N_SAMPLING_POINTS_KEY;
    extern const char* POINT_1_KEY;
    extern const char* POINT_2_KEY;
    extern const char* SURFACE_INTERSECTION_KEY;
    extern const char* FEATURE_LINE_ID_KEY;
    extern const char* FILE_NAME_KEY;
}

namespace fieldSamplingFunctionObjectKeys
{
    extern const char* SAMPLING_FILE_KEY;
}

namespace sceneKeys
{
    extern const char* SCENES_DICT_KEY;
    extern const char* BACKGROUND_COLORS_DICT_KEY;
    extern const char* WIDGETS_DICT_KEY;
    extern const char* CAMERAS_DICT_KEY;
    extern const char* ITEMS_DICT_KEY;
    extern const char* COLOR_LEGEND_OPTIONS_DICT_KEY;
    extern const char* COLOR_LEGENDS_DICT_KEY;
}

namespace backgroundKeys
{
    extern const char* BACKGROUND_1_KEY;
    extern const char* BACKGROUND_2_KEY;
}

namespace cameraKeys
{
    extern const char* FOCAL_POINT_KEY;
    extern const char* POSITION_KEY;
    extern const char* PARALLEL_PROJECTION_KEY;
    extern const char* UP_KEY;
    extern const char* PARALLEL_SCALE_KEY;
    extern const char* NAME_KEY;
    extern const char* COORDINATE_TYPE_KEY;
    extern const char* REFERENCE_FRAME_KEY;
}

namespace axisWidgetKeys
{
    extern const char* COLOUR_KEY;
    extern const char* VISIBLE_KEY;
}

namespace gridWidgetKeys
{
    extern const char* GRID_WIDGET_OPTIONS_DICT_KEY;
    extern const char* SHOW_AXES_GRID_KEY;
    extern const char* AXES_GRID_COLOR_KEY;
    extern const char* AXES_GRID_LABEL_COLOR_KEY;
    extern const char* SHOW_AXES_GRID_TITLE_KEY;
    extern const char* X_AXIS_TITLE_KEY;
    extern const char* Y_AXIS_TITLE_KEY;
    extern const char* Z_AXIS_TITLE_KEY;
    extern const char* MANUAL_FACE_SELECT_KEY;
    extern const char* SHOW_XYMIN_KEY;
    extern const char* SHOW_XYMAX_KEY;
    extern const char* SHOW_YZMIN_KEY;
    extern const char* SHOW_YZMAX_KEY;
    extern const char* SHOW_ZXMIN_KEY;
    extern const char* SHOW_ZXMAX_KEY;
    extern const char* USE_CUSTOM_BOUNDS_KEY;
    extern const char* X_BOUNDS_KEY;
    extern const char* Y_BOUNDS_KEY;
    extern const char* Z_BOUNDS_KEY;
    extern const char* USE_CUSTOM_LABELS_KEY;
    extern const char* USE_X_LABELS_KEY;
    extern const char* USE_Y_LABELS_KEY;
    extern const char* USE_Z_LABELS_KEY;
    extern const char* X_DISTANCE_LABEL_KEY;
    extern const char* Y_DISTANCE_LABEL_KEY;
    extern const char* Z_DISTANCE_LABEL_KEY;
}

namespace logoWidgetKeys
{
    extern const char* VISIBLE_KEY;
}

namespace timestepWidgetKeys
{
    extern const char* COLOUR_KEY;
    extern const char* VISIBLE_KEY;
    extern const char* FORMAT_KEY;
    extern const char* COORDINATES_KEY;
    extern const char* FONT_DICT_KEY;
}

namespace vectorWidgetKeys
{
    extern const char* VECTOR_WIDGET_OPTIONS_DICT_KEY;

    extern const char* VISIBLE_KEY;
    extern const char* VECTOR_REPRESENTATION_KEY;
    extern const char* VECTOR_3D_REP_OPTIONS_KEY;
    extern const char* VECTOR_TEXT_OPTIONS_KEY;

    extern const char* SHOW_VECTOR_3D_REP_KEY;
    extern const char* SHOW_VECTOR_RESULTANT_KEY;
    extern const char* SHOW_VECTOR_COMPONENTS_KEY;
    extern const char* VECTOR_WINDOW_SIZE_KEY;

    extern const char* SHOW_TITLE_KEY;
    extern const char* COORDINATES_KEY;

    extern const char* SHOW_VECTOR_TEXT_KEY;
    extern const char* LABEL_FORMAT_KEY;
    extern const char* FONT_SIZE_KEY;
    extern const char* SHOW_TEXT_BOX_KEY;

    extern const char* RUNTIME_INFO_DICT;
    extern const char* RESULTS_SUBDICT;
    extern const char* INTERFOAM_SUBDICT;
    extern const char* FLOWSOLVER_SUBDICT;
    extern const char* MAX_FRAME_ACCELERATION_COMPONENTS_KEY;
    extern const char* MIN_FRAME_ACCELERATION_COMPONENTS_KEY;
    extern const char* FRAME_ACCELERATION_KEY;
    extern const char* GRAVITY_ACCELERATION_KEY;

    extern const char* VECTOR_SUBDICT_KEY;
}

namespace referenceFramesWidgetKeys
{
extern const char *REFERENCE_FRAMES_WIDGET_LIST_KEY;
}

namespace activePatchKeys
{
    extern const char* ACTIVE_FIELD_KEY;
}

namespace meshObjectsKeys
{
    extern const char* TYPE_KEY;
    extern const char* REFERENCE_FRAME_ORIGIN_KEY;
    extern const char* REFERENCE_FRAME_E1_KEY;
    extern const char* REFERENCE_FRAME_E2_KEY;
    extern const char* REFERENCE_FRAME_E3_KEY;
    extern const char* REFERENCE_FRAME_TYPE_KEY;
    extern const char* MOTION_COORDINATE_FRAME_TYPE_KEY;
    extern const char* COORDINATE_SYSTEM_KEY;
    extern const char* GLOBAL_COORDINATE_SYSTEM_KEY;
    extern const char* INVALID_COORDINATE_SYSTEM_KEY;
}

namespace turbopropKeys
{
    extern const char* INLET_PATCHES_KEY;
    extern const char* OUTLET_PATCHES_KEY;
    extern const char* HUB_PATCHES_KEY;
    extern const char* SHROUD_PATCHES_KEY;
    extern const char* STREAM_ELEMENTS_KEY;
    extern const char* SPAN_ELEMENTS_KEY;
    extern const char* ORIGIN_KEY;
    extern const char* AXIS_KEY;
    extern const char* X_DIRECTION_KEY;
    extern const char* SPAN_VALUE_KEY;
    extern const char* STREAM_VALUE_KEY;
    extern const char* BLADE_TO_BLADE_KEY;
    extern const char* BLADE_PATCHES_KEY;
}

namespace groupKeys
{
    extern const char* SOURCES_KEY;
}

namespace thresholdKeys
{
    extern const char* FIELD_KEY;
    extern const char* ALL_POINTS_CRITERION_KEY;
    extern const char* MIN_THRESHOLD_KEY;
    extern const char* MAX_THRESHOLD_KEY;
    extern const char* INVERT_KEY;
}

namespace transformKeys
{
    extern const char* TRANSFORMATIONS_LIST_KEY;
    extern const char* TRANSFORMATION_TYPE_KEY;
    extern const char* INCLUDE_SOURCE_KEY;
    extern const char* COPY_COUNT_KEY;
    extern const char* OFFSET_DEGREES_KEY;
    extern const char* ORIGIN_KEY;
    extern const char* AXIS_KEY;
    extern const char* OFFSET_KEY;
    extern const char* SLICE_TYPE_DICT_KEY;
    extern const char* NORMAL_KEY;
    extern const char* DISTANCE_KEY;
    extern const char* ANGLE_DEGREES_KEY;
    extern const char* DIRECTION_KEY;
    extern const char* RATIO_KEY;
}

namespace importDataObjectKeys
{
    extern const char* TIME_STEP_KEY;
    extern const char* PATH_KEY;
    extern const char* FORMAT_KEY;
    extern const char* PVD_KEY;
    extern const char* VTK_KEY;
    extern const char* FIELD_SAMPLING_KEY;
    extern const char* ENSIGHT_KEY;
}

namespace objectKeys
{
    extern const char* SOURCE_KEY;
}

namespace tool3Dto2D
{
    extern const char* SOURCES_KEY;
    extern const char* ORIGIN_KEY;
    extern const char* NORMAL_KEY;
}

namespace surfaceSplitterKeys
{
    extern const char* ANGLE_KEY;
}

} // End namespace Foam
