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
    Foam::functionObjects::runTimeVis::ThresholdInfo

Description
    Information about the threshold objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/thresholdDataSetProvider.H"
#include "postDict/postDictKeys.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class ThresholdInfo : public ObjectInfo
{
protected:

    ThresholdInfo() = default;
    ThresholdObjectData data_;

public:
    ThresholdInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        ThresholdInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~ThresholdInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::THRESHOLD);

        data_.field = foamField(dict.lookup<word>(thresholdKeys::FIELD_KEY));
        data_.allPointsCriterion = dict.lookup<bool>(thresholdKeys::ALL_POINTS_CRITERION_KEY);
        data_.minThreshold = dict.lookup<scalar>(thresholdKeys::MIN_THRESHOLD_KEY);
        data_.maxThreshold = dict.lookup<scalar>(thresholdKeys::MAX_THRESHOLD_KEY);
        data_.invert = dict.lookup<bool>(thresholdKeys::INVERT_KEY);
    }

    [[nodiscard]] const ThresholdObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new ThresholdDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] ThresholdInfo* copy() const override {
        return new ThresholdInfo(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        return isTypeAndDataEqualTo<ThresholdInfo>(other, data_);
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        foamField thresholdField = foamField::withPointAssociation(data_.field);
        itemRequirements.addToRequiredFields(thresholdField);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
