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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ColourLegendsData

Description
    Contains all the colour legends data from the scene dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/widgets/colourLegendData.H"
#include "dictionaries/dictionaries.H"
#include "baseClasses/itemRequirements.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourLegendsData
{
    HashTable<ColourLegendData, word> colourLegends;

    // Construct from the scene dictionary
    explicit ColourLegendsData
    (
        const dictionary& sceneDict
    )
    {
        ColourLegendData defaultLegend = getDefaultColorLegend(sceneDict);

        const dictionary& legendsDict = GET_OPTIONAL_DICTIONARY(sceneDict, sceneKeys::COLOR_LEGENDS_DICT_KEY);

        colourLegends.clear();
        for(const std::string& legendName : legendsDict.toc())
        {
            foamField foamFieldName(legendName);
            ColourLegendData tempLegend
            (
                foamFieldName.lessAssociation(),
                legendsDict.subDict(legendName),
                defaultLegend,
                false
            );
            colourLegends.insert(foamFieldName.lessAssociation(), tempLegend);
        }
    }

    [[nodiscard]] static ColourLegendData getDefaultColorLegend(const dictionary &sceneDict)
    {
        if (sceneDict.found(sceneKeys::COLOR_LEGEND_OPTIONS_DICT_KEY))
        {
            const dictionary &defaultLegendsDict = sceneDict.subDict(sceneKeys::COLOR_LEGEND_OPTIONS_DICT_KEY);
            return {"default", defaultLegendsDict};
        }
        else
        {
            return {};
        }
    }

    [[nodiscard]] bool isAnyColorLegendVisible() const
    {
        return std::any_of(
            colourLegends.begin(),
            colourLegends.end(),
            [](const ColourLegendData &d) -> bool { return d.visible; }
        );
    }

    [[nodiscard]] ItemRequirements getItemRequirements() const
    {
        ItemRequirements requirements;
        for (const ColourLegendData& colourLegend : colourLegends)
        {
            if (colourLegend.visible)
            {
                requirements.addToRequiredFields(foamField::withCellAssociation(colourLegend.name));
            }
        }
        return requirements;
    }
};


} // End namespace
