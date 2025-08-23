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
    Foam::functionObjects::runTimeVis::FieldSamplingInfo

Description
    Information about the Field Sampling objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/fieldSamplingDataSetProvider.H"
#include "dataStructs/objects/fieldSamplingObjectData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class FieldSamplingInfo : public ObjectInfo
{
protected:

    FieldSamplingInfo() = default;
    FieldSamplingObjectData data_ {};

public:
    FieldSamplingInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict,
        const string& caseFolder
    )
    :
        ObjectInfo(dict)
    {
        FieldSamplingInfo::updateObjectData(dict);
        data_.caseFolder = caseFolder;
        readVisualisationDict(objectVisDict);
    };
    ~FieldSamplingInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::FIELD_SAMPLING);
        data_.readDict(dict);
    }

    [[nodiscard]] const FieldSamplingObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new FieldSamplingDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] FieldSamplingInfo* copy() const override {
        return new FieldSamplingInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        data_.computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        const auto* casted = dynamic_cast<const FieldSamplingInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        if (isVisible() && visualisation_.colourField.isCellAssociation())
        {
            foamField scalarField = foamField::withPointAssociation(visualisation_.colourField);
            itemRequirements.addToRequiredFields(scalarField);
        }
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
