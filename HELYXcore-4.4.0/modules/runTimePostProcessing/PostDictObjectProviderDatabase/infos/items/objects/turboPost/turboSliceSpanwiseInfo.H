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
    Foam::functionObjects::runTimeVis::TurboSliceSpanwiseInfo

Description
    Information about the Spanwise Slice Post-processing objects,
    which shows a slice with a constant stream.
    Contains the data from the dictionary.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/turboPost/turboSliceSpanwiseDataSetProvider.H"
#include "dataStructs/objects/turboPost/turboSliceSpanwiseObjectData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class TurboSliceSpanwiseInfo : public ObjectInfo
{
protected:

    TurboSliceSpanwiseInfo() = default;
    TurboSliceSpanwiseObjectData data_ {};

public:
    TurboSliceSpanwiseInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        TurboSliceSpanwiseInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~TurboSliceSpanwiseInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::TURBO_SLICE_SPANWISE);
        data_.readDict(dict);
    }

    [[nodiscard]] const TurboSliceSpanwiseObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new TurboSliceSpanwiseDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] TurboSliceSpanwiseInfo* copy() const override {
        return new TurboSliceSpanwiseInfo(*this);
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
        const auto* casted = dynamic_cast<const TurboSliceSpanwiseInfo*>(other);
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
