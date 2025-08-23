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
    Foam::functionObjects::runTimeVis::ColourLookupTableInfo

Description
    Information about the lookup tables
    Contains the data from the dictionary

SourceFiles
    <colourLookupTableInfo.C>

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/scene/colourLookupTableData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ColourLookupTablesInfo
{
private:

    HashTable<ColourLookupTableData, word> colourLutData_;
    ColourLookupTableData defaultLut_;

    [[nodiscard]] ColourLookupTableData getColorLookupTableData(const word& fieldLA) const;

public:
    explicit ColourLookupTablesInfo
    (
            const dictionary& dict
    );

    [[nodiscard]] ColourLookupTableData getColorLookupTableData
    (
        const foamField& field,
        const ColourLookupTablesInfo& baseFields
    ) const;

    [[nodiscard]] List<word> getTableOfContents() const {return colourLutData_.toc();};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam

// ************************************************************************* //
