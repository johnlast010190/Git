/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::Tool3Dto2DInfo

Description
    Information about the 3Dto2D objects

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "itemDataSetProviders/testObjects/tool3Dto2DDataSetProvider.H"
#include "dataStructs/objects/tool3Dto2DObjectData.H"

#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class Tool3Dto2DInfo : public ItemInfo
{
protected:

    Tool3Dto2DInfo() : ItemInfo(Visualisation::solidColorDefaults()) {};
    Tool3Dto2DObjectData data_ {};

public:
    Tool3Dto2DInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ItemInfo(dict, ItemType::Value::TOOL_3D_TO_2D, Visualisation::solidColorDefaults())
    {
        Tool3Dto2DInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~Tool3Dto2DInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::TOOL_3D_TO_2D);
        data_.readDict(dict);
    }

    [[nodiscard]] const Tool3Dto2DObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new Tool3Dto2DDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] Tool3Dto2DInfo* copy() const override {
        return new Tool3Dto2DInfo(*this);
    }

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        return data_.sourceSurfaces;
    };

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        data_.computeAndAddToHash(hash);
        return hash;
    }

    [[nodiscard]] bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const Tool3Dto2DInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override
    {
        return false;
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
