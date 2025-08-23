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
    Foam::functionObjects::runTimeVis::BoundariesGroupInfo

Description
    Created to combine all the external/internal boundaries of a region

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

class BoundariesGroupInfo : public ItemInfo
{
protected:
    BoundariesGroupInfo() : ItemInfo(Visualisation::solidColorDefaults()) {};
    GroupData data_;

public:
    BoundariesGroupInfo
    (
        const std::string& name,
        const std::string& region,
        const std::vector<Id>& patchIds
    )
    : ItemInfo(Visualisation::solidColorDefaults())
    {
        id_.name = name;
        id_.region = region;
        id_.type = ItemType(ItemType::Value::BOUNDARY_GROUP);

        for (const Id& memberId : patchIds)
        {
            data_.sources.push_back(memberId);
        }
    };

    ~BoundariesGroupInfo() override = default;

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        return data_.sources;
    };

    const GroupData* getData() { return &data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new GroupDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] BoundariesGroupInfo* copy() const override {
        return new BoundariesGroupInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return 0;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const BoundariesGroupInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
