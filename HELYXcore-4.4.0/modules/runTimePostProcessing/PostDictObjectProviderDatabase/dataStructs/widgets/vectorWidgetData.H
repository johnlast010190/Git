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
    Foam::functionObjects::runTimeVis::VectorWidgetData

Description
    Information about the vector widget
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point2D/point2D.H"
#include "types/vectorDataType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct VectorWidget3dData
{
    bool showVector3D = false;
    bool showResultant = true;
    bool showComponents = false;
    scalar size = 0.25;
    bool showTitle = false;
    point2D coordinates = {0.05, 0.8};

    VectorWidget3dData() = default;

    explicit VectorWidget3dData(const dictionary& vectorDict)
    {
        if (vectorDict.found(vectorWidgetKeys::VECTOR_3D_REP_OPTIONS_KEY))
        {
            const dictionary &vector3dDict = vectorDict.subDict(vectorWidgetKeys::VECTOR_3D_REP_OPTIONS_KEY);
            showVector3D = vector3dDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_VECTOR_3D_REP_KEY, false);
            showResultant = vector3dDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_VECTOR_RESULTANT_KEY, true);
            showComponents = vector3dDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_VECTOR_COMPONENTS_KEY, false);
            size = vector3dDict.lookupOrDefault<scalar>(vectorWidgetKeys::VECTOR_WINDOW_SIZE_KEY, 1) / 4;
            showTitle = vector3dDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_TITLE_KEY, false);
            coordinates = vector3dDict.lookupOrDefault<point2D>(vectorWidgetKeys::COORDINATES_KEY, {0.05, 0.8});
        }
    }
};

struct VectorWidgetTextData
{
    bool showVectorTextRep = false;
    string labelFormat = "%.2f";
    label fontSize = 12;
    bool showTextBox = false;
    bool showTitle = false;
    point2D coordinates = {0.125, 0.8};

    VectorWidgetTextData() = default;

    explicit VectorWidgetTextData(const dictionary& vectorDict)
    {
        if (vectorDict.found(vectorWidgetKeys::VECTOR_TEXT_OPTIONS_KEY))
        {
            const dictionary &textDict = vectorDict.subDict(vectorWidgetKeys::VECTOR_TEXT_OPTIONS_KEY);
            showVectorTextRep = textDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_VECTOR_TEXT_KEY, false);
            labelFormat = textDict.lookupOrDefault<string>(vectorWidgetKeys::LABEL_FORMAT_KEY, "%.2f");
            fontSize = textDict.lookupOrDefault<label>(vectorWidgetKeys::FONT_SIZE_KEY, 12);
            showTextBox = textDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_TEXT_BOX_KEY, false);
            showTitle = textDict.lookupOrDefault<bool>(vectorWidgetKeys::SHOW_TITLE_KEY, false);
            coordinates = textDict.lookupOrDefault<point2D>(vectorWidgetKeys::COORDINATES_KEY, {0.125, 0.8});
        }
    }
};

struct VectorWidgetData
{
protected:
    const dictionary& vectorSubDict_;

public:
    bool visible = false;
    VectorDataType dataSource;
    VectorWidget3dData vector3dData;
    VectorWidgetTextData vectorTextData;

    explicit VectorWidgetData(const dictionary& coloursDict)
    :
        vectorSubDict_(coloursDict.optionalSubDict(vectorWidgetKeys::VECTOR_WIDGET_OPTIONS_DICT_KEY)),
        vector3dData(vectorSubDict_),
        vectorTextData(vectorSubDict_)
    {
        visible = vectorSubDict_.lookupOrDefault<bool>(vectorWidgetKeys::VISIBLE_KEY, false);
        dataSource = vectorSubDict_.lookupOrDefault(vectorWidgetKeys::VECTOR_REPRESENTATION_KEY, VectorDataType());
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
