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
    Foam::functionObjects::runTimeVis::BoundarySurfaceInfo

Description
    Information about a user defined block mesh boundary surface
    Contains the data from the blockMeshDict dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "itemDataSetProviders/geometry/boundarySurfaceDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class BoundarySurfaceInfo : public ItemInfo
{
protected:
    BoundarySurfaceData data_;

public:
    explicit BoundarySurfaceInfo
    (
        const BlockMeshPatch& patch,
        const BlockMeshBlock& block,
        const List<point>& vertices
    )
    :
        ItemInfo(Visualisation::indexDefaults()),
        data_(patch, block, vertices)
    {
        id_.name = patch.name;
        id_.type = ItemType::Value::SURFACE;
        id_.region = regionKeys::DEFAULT_REGION_KEY;
    };

    ~BoundarySurfaceInfo() override = default;

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new BoundarySurfaceDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] BoundarySurfaceInfo* copy() const override {
        return new BoundarySurfaceInfo(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        auto casted = dynamic_cast<const BoundarySurfaceInfo*>(other);
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

} // End namespace

// ************************************************************************* //
