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
    Foam::functionObjects::runTimeVis::ClipInfo

Description
    Information about the clip objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "infos/items/objects/cutting/cuttingInfo.H"
#include "itemDataSetProviders/objects/clipDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class ClipInfo : public ObjectInfo
{
protected:

    ClipObjectData data_ {};
    CuttingInfo cuttingInfo_ {};
    string caseFolder_;
    const ReferenceFrames& referenceFrames_;

public:
    ClipInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict,
        string caseFolder,
        const ReferenceFrames& referenceFrames
    )
    :
        ObjectInfo(dict),
        referenceFrames_(referenceFrames)
    {
        caseFolder_ = std::move(caseFolder);
        ClipInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~ClipInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::CLIP);
        const dictionary& sliceTypeDict = dict.subDict(clipObjectKeys::SLICE_TYPE_DICT_KEY);
        cuttingInfo_.readDict(sliceTypeDict, caseFolder_);
        data_.readDict(dict);
    }

    [[nodiscard]] const ClipObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new ClipDataSetProvider(id_.name, data_, cuttingInfo_.createCuttingSurfaceProvider(), referenceFrames_);
    };

    [[nodiscard]] ClipInfo* copy() const override {
        return new ClipInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = data_.computeAndReturnHash();
        cuttingInfo_.computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const ClipInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_ &&
               cuttingInfo_ == casted->cuttingInfo_;
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        cuttingInfo_.addItemRequirements(itemRequirements);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
