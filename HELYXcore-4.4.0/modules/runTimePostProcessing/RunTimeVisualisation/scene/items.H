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
    Foam::functionObjects::runTimeVis::Items

Description
    A class that holds a collection of items

SourceFiles
    items.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "scene/item.H"
#include "colourLookupTable/colourLookupTableProvider.H"
#include "rendering/renderManager.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class PostDictObjectProviderDatabase;

/*---------------------------------------------------------------------------*\
                          Class Items Declaration
\*---------------------------------------------------------------------------*/

class Items
{
private:
    std::vector<Item> items_;

    const PostDictObjectProviderDatabase &database_;
    const std::string sceneName_;

public:

    // Constructors

    Items(
        const PostDictObjectProviderDatabase &database,
        const HashTable<const ItemInfo *, Id, IdHasher> &itemInfos,
        std::string sceneName
    );

    // Member Functions

    void updateActors
        (
            label currentTimeLabel,
            scalar currentTime,
            ColourLookupTableProvider &provider,
            const ColourMaps &colourMaps
        );
    vtkSmartPointer<vtkDataSet> createDataSetWithAllData(label currentTimeLabel, scalar currentTime);

    void registerDataToCompositer(Compositer *compositer);
    void redistributeDataWithCompositerIfNecessary(Compositer *compositer);

    void addToRenderer(renderManager &renderer);

    [[nodiscard]] bool isAnyVisible() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
