/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::InvalidItemInfo

Description
    Information about items that are requested but have some issue

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "itemDataSetProviders/missingItem/missingItemDataSetProvider.H"
#include "postDict/postDictKeys.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class InvalidItemInfo : public ItemInfo
{
public:
    explicit InvalidItemInfo(const Id& id)
    :
        ItemInfo(Visualisation::solidColorDefaults())
    {
        id_.type = id.type;
        id_.name = id.name;
        id_.region = id.region;
    };

    InvalidItemInfo(const string& name, const ItemType& type)
        :
        ItemInfo(Visualisation::solidColorDefaults())
    {
        id_.type = type;
        id_.name = name;
        id_.region = regionKeys::DEFAULT_REGION_KEY;
    };

    ~InvalidItemInfo() override = default;

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        return {};
    };

    [[nodiscard]] bool isValid() const override {return false;};
    [[nodiscard]] word getInvalidReason() const override = 0;

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new MissingItemDataSetProvider(id_.name);
    };

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return 0;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const InvalidItemInfo*>(other);
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

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
