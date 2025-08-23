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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::StreamlinesInfo

Description
    Information about the streamline objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "infos/items/objects/streamlineSources/streamlineSourceInfo.H"
#include "itemDataSetProviders/objects/streamlinesDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class StreamlinesInfo : public ObjectInfo
{
protected:

    StreamlinesInfo() = default;
    StreamlinesObjectData data_;
    StreamlineSourceInfo pointSourceInfo_;

public:
    StreamlinesInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        StreamlinesInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~StreamlinesInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::STREAMLINES);
        pointSourceInfo_.readDict(dict);
        data_.readDict(dict);
    }

    [[nodiscard]] const StreamlinesObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new StreamlinesDataSetProvider(
            id_.name,
            data_,
            pointSourceInfo_.createStreamlineSourceProvider());
    };

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        std::vector<Id> sources;
        sources.push_back(source_);
        pointSourceInfo_.addToSources(sources);
        return sources;
    }

    [[nodiscard]] StreamlinesInfo* copy() const override {
        return new StreamlinesInfo(*this);
    }


protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        data_.computeAndAddToHash(hash);
        pointSourceInfo_.computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const StreamlinesInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_ &&
               pointSourceInfo_ == casted->pointSourceInfo_;
    };
    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        foamField vectorField = foamField::withPointAssociation(data_.vectorField);
        itemRequirements.addToRequiredFields(vectorField);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
