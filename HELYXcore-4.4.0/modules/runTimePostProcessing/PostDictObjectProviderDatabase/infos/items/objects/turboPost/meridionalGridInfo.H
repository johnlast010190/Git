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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MeridionalGridInfo

Description
    Information about the Turbo Background Mesh Post-processing objects,
    which is the basis for the other 3. Contains the data from the dictionary.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/turboPost/meridionalGridDataSetProvider.H"
#include "dataStructs/objects/turboPost/meridionalGridObjectData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class MeridionalGridInfo : public ObjectInfo
{
protected:

    MeridionalGridInfo() = default;
    MeridionalGridObjectData data_ {};

public:
    MeridionalGridInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        MeridionalGridInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~MeridionalGridInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::MERIDIONAL_GRID);
        data_.readDict(dict);
    }

    [[nodiscard]] const MeridionalGridObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new MeridionalGridDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] MeridionalGridInfo* copy() const override {
        return new MeridionalGridInfo(*this);
    }

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        std::vector<Id> sources;
        sources.push_back(source_);
        addSourcesFrom(data_.inletPatches, sources);
        addSourcesFrom(data_.outletPatches, sources);
        addSourcesFrom(data_.hubPatches, sources);
        addSourcesFrom(data_.shroudPatches, sources);
        return sources;
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
        const auto* casted = dynamic_cast<const MeridionalGridInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

private:
    void addSourcesFrom(const List<word>& patches, std::vector<Id>& sources) const
    {
        for (const word& inletPatch : patches)
        {
            sources.emplace_back(inletPatch, getId().region, ItemType::Value::PATCH);
        }
    }

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
