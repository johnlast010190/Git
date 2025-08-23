/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::GridWidgetData

Description
    Contains the dict data for the grid widget

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "primitives/Vector/Vector.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct GridWidgetData
{
    Vector<scalar> axesGridColor;
    Vector<scalar> axesGridLabelColor;
    bool showAxesGrid;
    bool showAxesGridTitle;
    string XAxisTitle;
    string YAxisTitle;
    string ZAxisTitle;
    bool manualFaceSelect;
    bool showXYMin;
    bool showXYMax;
    bool showYZMin;
    bool showYZMax;
    bool showZXMin;
    bool showZXMax;
    bool useCustomBounds;
    Pair<scalar> XBounds {};
    Pair<scalar> YBounds {};
    Pair<scalar> ZBounds {};
    bool useCustomLabels;
    scalar XDistanceLabel;
    scalar YDistanceLabel;
    scalar ZDistanceLabel;
    bool useXDistanceLabel;
    bool useYDistanceLabel;
    bool useZDistanceLabel;

    explicit GridWidgetData(const dictionary& coloursDict)
    {
        static const Vector<scalar> WHITE = {1,1,1};
        static const Vector<scalar> GRAY = {0.5,0.5,0.5};

        const dictionary& widgetDict = coloursDict.optionalSubDict(gridWidgetKeys::GRID_WIDGET_OPTIONS_DICT_KEY);
        showAxesGrid  = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_AXES_GRID_KEY, false);
        axesGridColor = widgetDict.lookupOrDefault<Vector<scalar>>(gridWidgetKeys::AXES_GRID_COLOR_KEY, WHITE);
        axesGridLabelColor = widgetDict.lookupOrDefault<Vector<scalar>>(gridWidgetKeys::AXES_GRID_LABEL_COLOR_KEY, GRAY);
        showAxesGridTitle = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_AXES_GRID_TITLE_KEY, false);
        XAxisTitle = widgetDict.lookupOrDefault<string>(gridWidgetKeys::X_AXIS_TITLE_KEY, "X Axis");
        YAxisTitle = widgetDict.lookupOrDefault<string>(gridWidgetKeys::Y_AXIS_TITLE_KEY, "Y Axis");
        ZAxisTitle = widgetDict.lookupOrDefault<string>(gridWidgetKeys::Z_AXIS_TITLE_KEY, "Z Axis");
        manualFaceSelect = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::MANUAL_FACE_SELECT_KEY, false);
        showXYMin = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_XYMIN_KEY, false);
        showXYMax = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_XYMAX_KEY, false);
        showYZMin = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_YZMIN_KEY, false);
        showYZMax = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_YZMAX_KEY, false);
        showZXMin = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_ZXMIN_KEY, false);
        showZXMax = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::SHOW_ZXMAX_KEY, false);
        useCustomBounds = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::USE_CUSTOM_BOUNDS_KEY, false);
        XBounds = widgetDict.lookupOrDefault<Pair<scalar>>(gridWidgetKeys::X_BOUNDS_KEY, {0,1});
        YBounds = widgetDict.lookupOrDefault<Pair<scalar>>(gridWidgetKeys::Y_BOUNDS_KEY, {0,1});
        ZBounds = widgetDict.lookupOrDefault<Pair<scalar>>(gridWidgetKeys::Z_BOUNDS_KEY, {0,1});
        useCustomLabels = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::USE_CUSTOM_LABELS_KEY, false);
        useXDistanceLabel = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::USE_X_LABELS_KEY, false);
        useYDistanceLabel = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::USE_Y_LABELS_KEY, false);
        useZDistanceLabel = widgetDict.lookupOrDefault<bool>(gridWidgetKeys::USE_Z_LABELS_KEY, false);
        XDistanceLabel = widgetDict.lookupOrDefault<scalar>(gridWidgetKeys::X_DISTANCE_LABEL_KEY, 0);
        YDistanceLabel = widgetDict.lookupOrDefault<scalar>(gridWidgetKeys::Y_DISTANCE_LABEL_KEY, 0);
        ZDistanceLabel = widgetDict.lookupOrDefault<scalar>(gridWidgetKeys::Z_DISTANCE_LABEL_KEY, 0);
    }
};

} // End namespace

// ************************************************************************* //
