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
    Foam::functionObjects::runTimeVis::ColourLegendData

Description
    Contains the colour legend data from the dictionary

SourceFiles
    colorLegendData.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/word/word.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/Vector/Vector.H"
#include "primitives/ranges/MinMax/MinMax.H"

#include "baseClasses/foamField.H"
#include "postDict/postDictKeys.H"
#include "types/location.H"
#include "infos/fontInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourLegendData
{
    foamField name;
    string title;
    bool visible;
    string labelFormat;
    label numberOfLabels;
    Location location;
    Tuple2<scalar, scalar> coordinates;
    label thickness;
    scalar length;
    bool vertical;
    bool autoLabels;
    bool showTicks;

    FontInfo font;

    // Construct from a dictionary
    ColourLegendData
    (
        const word& name,
        const dictionary& legendDict
    );
    // Construct from a dictionary, using the given default values
    ColourLegendData
    (
        const word& name,
        const dictionary& legendDict,
        const ColourLegendData& defaultColorLegendData,
        bool allowNoTitle
    );

    // Just default
    ColourLegendData();
};


} // End namespace
