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

\*---------------------------------------------------------------------------*/

#include <utility>

#include "itemInfoIndex.H"
#include "itemStorage.H"

namespace Foam::functionObjects::runTimeVis
{

ItemInfoIndex::ItemInfoIndex(const ItemInfo* info, const ItemStorage* storage, std::string  sceneName)
    : itemInfo(info),
      storage(storage),
      sceneName(std::move(sceneName))
{
    hash = getHashForItemInfoAndSources(itemInfo);
}

bool ItemInfoIndex::hasEqualSource(const ItemInfoIndex &otherSource) const {
    bool found = false;
    for (const Id &id : itemInfo->getSources()) {
        ItemInfoIndex source(storage->getItemInfo(id, sceneName), storage, sceneName);
        if (source == otherSource) {
            found = true;
            break;
        }
    }
    return found;
}

bool ItemInfoIndex::operator==(const ItemInfoIndex& other) const
{
    if (!itemInfo->isDataEqualTo(other.itemInfo))
    {
        return false;
    }
    if (itemInfo->getSources().size() != other.itemInfo->getSources().size())
    {
        return false;
    }
    std::vector<Id> sourceIds = itemInfo->getSources();
    return std::all_of(sourceIds.begin(), sourceIds.end(), [this, &other](const Id& id)
        {
                ItemInfoIndex source(storage->getItemInfo(id, sceneName), storage, sceneName);
                return other.hasEqualSource(source);
        });
}

bool ItemInfoIndex::operator!=(const ItemInfoIndex& other) const
{
    return !(operator==(other));
}

size_t ItemInfoIndex::getHashForItemInfoAndSources(const ItemInfo* info) const
{
    size_t calcHash = info->computeAndReturnItemInfoHash();
    size_t res;
    for (const Id &id : info->getSources()) {
        res = getHashForItemInfoAndSources(storage->getItemInfo(id, sceneName));
        hasher::hash_add_unordered(calcHash, res);
    }
    return calcHash;
}

// ************************************************************************* //
} // End namespace
