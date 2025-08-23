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
    Foam::functionObjects::runTimeVis::ColourLookupTableData

Description
    Contains the colour map configuration data from the dictionary for a given scene

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/ranges/MinMax/MinMax.H"
#include "postDict/postDictKeys.H"
#include "baseClasses/foamField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourLookupTableData
{
    foamField fieldName {string("")};
    word colourMapName = "ENGYS_DEFAULT";
    label resolution = 256;
    bool inverted = false;
    bool automaticRange = true;
    scalarMinMax fixedRange = Foam::scalarMinMax(42, -42);

    // Construct from a dictionary
    ColourLookupTableData
    (
        const foamField& name,
        const dictionary& colourLegendDict
    )
    :
        ColourLookupTableData(name, colourLegendDict, ColourLookupTableData())
    {};

    // Construct from a dictionary with default values
    ColourLookupTableData
    (
        const foamField& name,
        const dictionary& colourLegendDict,
        const ColourLookupTableData& defaultMap
    )
    :
        fieldName(name)
    {
        readDict(colourLegendDict, defaultMap);
    }

    ColourLookupTableData() = default;


private:

    void readDict
    (
        const dictionary& colourLegendDict,
        const ColourLookupTableData& defaultMap
    )
    {
        if (colourLegendDict.found(colourLookupTableKeys::COLOUR_LOOKUP_DICT_KEY))
        {
            const dictionary& colourMapDict = colourLegendDict.subDict(colourLookupTableKeys::COLOUR_LOOKUP_DICT_KEY);

            colourMapName = colourMapDict.lookupOrDefault<word>(colourLookupTableKeys::TYPE_KEY, defaultMap.colourMapName);
            resolution = colourMapDict.lookupOrDefault<label>(colourLookupTableKeys::RESOLUTION_KEY, defaultMap.resolution);
            inverted = colourMapDict.lookupOrDefault<bool>(colourLookupTableKeys::INVERTED_KEY, defaultMap.inverted);
            automaticRange = colourMapDict.lookupOrDefault<bool>(colourLookupTableKeys::AUTOMATIC_RANGE_KEY, defaultMap.automaticRange);
            fixedRange = colourMapDict.lookupOrDefault<scalarMinMax>(colourLookupTableKeys::RANGE_KEY, defaultMap.fixedRange);
        }
        else
        {
            colourMapName = defaultMap.colourMapName;
            resolution = defaultMap.resolution;
            inverted = defaultMap.inverted;
            automaticRange = defaultMap.automaticRange;
            fixedRange = defaultMap.fixedRange;
        }
    }
};


} // End namespace
