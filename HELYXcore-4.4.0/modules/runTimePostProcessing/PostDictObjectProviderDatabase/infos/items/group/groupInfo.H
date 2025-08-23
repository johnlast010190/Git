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
    Foam::functionObjects::runTimeVis::GroupInfo

Description
    Information about the group objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "itemDataSetProviders/group/groupDataSetProvider.H"
#include "postDict/postDictKeys.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class GroupInfo : public ItemInfo
{
protected:
    GroupInfo() : ItemInfo(Visualisation::solidColorDefaults()) {};
    GroupData data_;

public:
    explicit GroupInfo
    (
        const dictionary& itemDict
    )
    :
        ItemInfo(itemDict, ItemType::Value::GROUP, Visualisation::solidColorDefaults())
    {
        GroupInfo::updateObjectData(itemDict);
    };

    void updateObjectData(const dictionary& itemDict) override
    {
        if (data_.sources.empty() || itemDict.found(groupKeys::SOURCES_KEY))
        {
            data_.sources.clear();
            auto sources = itemDict.lookup<List<dictionary>>(groupKeys::SOURCES_KEY);
            for (const dictionary &d: sources)
            {
                data_.sources.emplace_back(d);
            }
        }
    }

    ~GroupInfo() override = default;

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        return data_.sources;
    };

    const GroupData* getData() { return &data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new GroupDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] GroupInfo* copy() const override {
        return new GroupInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return 0;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const GroupInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
