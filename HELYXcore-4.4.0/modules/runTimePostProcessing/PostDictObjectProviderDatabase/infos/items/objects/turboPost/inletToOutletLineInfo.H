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
    Foam::functionObjects::runTimeVis::InletToOutletLineInfo

Description
    Information about the Inlet-to-Outlet Post-processing objects,
    which shows a line along a constant span.
    Contains the data from the dictionary.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/turboPost/inletToOutletLineDataSetProvider.H"
#include "dataStructs/objects/turboPost/inletToOutletLineObjectData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class InletToOutletLineInfo : public ObjectInfo
{
protected:

    InletToOutletLineInfo() = default;
    InletToOutletLineObjectData data_ {};

public:
    InletToOutletLineInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        InletToOutletLineInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~InletToOutletLineInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::INLET_TO_OUTLET_LINE);
        data_.readDict(dict);
    }

    [[nodiscard]] const InletToOutletLineObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new InletToOutletLineDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] InletToOutletLineInfo* copy() const override {
        return new InletToOutletLineInfo(*this);
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
        const auto* casted = dynamic_cast<const InletToOutletLineInfo*>(other);
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
