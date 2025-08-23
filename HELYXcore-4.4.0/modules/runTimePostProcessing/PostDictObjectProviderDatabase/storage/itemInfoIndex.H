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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ItemInfoIndex

Description
    Wrapper for an ItemInfo that allows them to be used as a key in a HashTable

SourceFiles
    itemInfoIndex.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
class ItemStorage;

struct ItemInfoIndex
{
    const ItemInfo* itemInfo;
    const ItemStorage* storage;
    std::string sceneName;
    size_t hash;

    ItemInfoIndex(const ItemInfo* info, const ItemStorage* storage, std::string  sceneName);

    ItemInfoIndex() {itemInfo = nullptr; storage = nullptr; hash = 0;};

    bool operator==(const ItemInfoIndex& other) const;

    bool operator!=(const ItemInfoIndex& other) const;

    friend Ostream& operator<<(Ostream& os, const ItemInfoIndex& info)
    {
        os << info.itemInfo->getId();
        return os;
    };

private:
    [[nodiscard]] bool hasEqualSource(const ItemInfoIndex &source) const;

    size_t getHashForItemInfoAndSources(const ItemInfo* info) const;
};

struct ItemInfoIndexHasher
{
    std::size_t operator()(const ItemInfoIndex& key) const {return key.hash;};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
