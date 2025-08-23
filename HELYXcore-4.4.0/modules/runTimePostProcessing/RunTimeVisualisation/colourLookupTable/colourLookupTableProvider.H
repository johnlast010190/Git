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

Description
    Class that generates VTK colour lookup tables as required for the used fields

Class
    Foam::functionObjects::ColourLookupTableProvider

\*---------------------------------------------------------------------------*/

#pragma once

#include "colourLookupTable/colourLookupTable.H"
#include "infos/colourLookupTablesInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourMaps;

/*---------------------------------------------------------------------------*\
               Class ColourLookupTableProvider Declaration
\*---------------------------------------------------------------------------*/

class ColourLookupTableProvider
{
protected:

    const ColourLookupTablesInfo sceneLutInfo_;
    const ColourLookupTablesInfo& baseLutInfo_;

    HashTable<ColourLookupTable, word> colourLookupTables;

public:

    // Constructors

    ColourLookupTableProvider
        (
            ColourLookupTablesInfo sceneLutInfo,
            const ColourLookupTablesInfo& baseLutInfo
        );


    // Member Functions

    const ColourLookupTable *updateAndReturnColorLookupTable
        (
            const foamField &field,
            const ColourMaps &colourMaps,
            const scalarMinMax &domainMinMax
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
