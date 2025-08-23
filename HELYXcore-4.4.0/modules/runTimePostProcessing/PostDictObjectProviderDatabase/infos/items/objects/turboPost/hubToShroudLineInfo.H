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
    (c) 2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::HubToShroudLineInfo

Description
    Information about the Hub-to-Shroud Post-processing objects,
    which shows a line along a constant stream.
    Contains the data from the dictionary.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/turboPost/hubToShroudLineDataSetProvider.H"
#include "dataStructs/objects/turboPost/hubToShroudLineObjectData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class HubToShroudLineInfo : public ObjectInfo
{
protected:

    HubToShroudLineInfo() = default;
    HubToShroudLineObjectData data_ {};

public:
    HubToShroudLineInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        HubToShroudLineInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~HubToShroudLineInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::HUB_TO_SHROUD_LINE);
        data_.readDict(dict);
    }

    [[nodiscard]] const HubToShroudLineObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new HubToShroudLineDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] HubToShroudLineInfo* copy() const override {
        return new HubToShroudLineInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        data_.computeAndAddToHash(hash);
        return hash;
    }

    [[nodiscard]] bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const HubToShroudLineInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
