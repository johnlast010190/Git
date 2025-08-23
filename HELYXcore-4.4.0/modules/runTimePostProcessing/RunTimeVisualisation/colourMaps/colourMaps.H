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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ColorMaps

Description
    A struct with all available color maps, defined in the code

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "colourMaps/colourMap.H"
#include "colourMaps/availableColourmaps.H"
#include "db/IOstreams/StringStreams/IStringStream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourMaps {

    // Stores definitions of colourmaps (e.g. rgb values for Viridis,
    // RedToBlue, etc...)
    HashTable<ColourMap, word> colourMaps_;

    ColourMaps()
    {
        // availableColourMapsCharArray defined in the availableColourMaps header,
        // which is written by CMake at generation-time
        IStringStream temp(availableColourMapsCharArray);
        dictionary availableColourMapsDict(temp);
        availableColourMapsDict =
            &availableColourMapsDict.subDict("availableColourMaps");

        for (const word& key : availableColourMapsDict.toc())
        {
            ColourMap cmInfo(availableColourMapsDict.subDict(key));
            colourMaps_.set(cmInfo.name, cmInfo);
        }
    };

    [[nodiscard]] const ColourMap* getColourMap(const word& colormap) const
    {
        return &colourMaps_[colormap];
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
