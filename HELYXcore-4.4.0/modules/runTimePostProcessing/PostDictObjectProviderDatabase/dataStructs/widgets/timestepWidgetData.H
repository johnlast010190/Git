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
    Foam::functionObjects::runTimeVis::TimestepWidgetData

Description
    Information about the timestep widget
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/point.H"
#include "meshes/primitiveShapes/point2D/point2D.H"
#include "infos/fontInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct TimestepWidgetData
{
    bool visible = false;
    string format;
    point2D coordinates;
    FontInfo font;

    TimestepWidgetData() = default;

    explicit TimestepWidgetData(const dictionary& coloursDict)
    {
        visible = coloursDict.lookupOrDefault<bool>(timestepWidgetKeys::VISIBLE_KEY, false);
        format = coloursDict.lookupOrDefault<string>(timestepWidgetKeys::FORMAT_KEY, "%6.2f [s]");
        coordinates = coloursDict.lookupOrDefault<point2D>(timestepWidgetKeys::COORDINATES_KEY, point2D(0.95,0.95));

        FontInfo defaultFont;
        defaultFont.size = 20;
        if (coloursDict.found(timestepWidgetKeys::FONT_DICT_KEY))
        {
            const dictionary& labelFontDict = &coloursDict.subDict(timestepWidgetKeys::FONT_DICT_KEY);
            font.readDict(labelFontDict, defaultFont);
        }
        else
        {
            font = defaultFont;
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
