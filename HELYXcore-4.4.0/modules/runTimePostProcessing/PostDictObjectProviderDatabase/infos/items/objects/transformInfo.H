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
    Foam::functionObjects::runTimeVis::TransformInfo

Description
    Information about the transform objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/transformDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class TransformInfo : public ObjectInfo
{
protected:

    TransformInfo() = default;
    TransformObjectData data_;

public:
    TransformInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        TransformInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~TransformInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::TRANSFORM);
        data_.transformations.clear();
        auto transformationDicts = dict.lookup<List<dictionary>>(transformKeys::TRANSFORMATIONS_LIST_KEY);
        for (dictionary d : transformationDicts)
        {
            SingleTransformationObjectData singleTransformation;

            singleTransformation.type = d.lookup<TransformationType>(transformKeys::TRANSFORMATION_TYPE_KEY);
            singleTransformation.includeSource = d.lookupOrDefault(transformKeys::INCLUDE_SOURCE_KEY, false);
            switch (singleTransformation.type.getValue()) {
            case TransformationType::CYLINDRICAL_ARRAY:
                singleTransformation.transformCount = d.lookup<label>(transformKeys::COPY_COUNT_KEY);
                singleTransformation.transformOffset = d.lookup<scalar>(transformKeys::OFFSET_DEGREES_KEY);
                singleTransformation.transformPoint = d.lookup<point>(transformKeys::ORIGIN_KEY);
                singleTransformation.transformVector = d.lookup<vector>(transformKeys::AXIS_KEY);
                break;
            case TransformationType::MIRROR: {
                singleTransformation.transformCount = 0;
                singleTransformation.transformOffset = d.lookup<scalar>(transformKeys::OFFSET_KEY);
                dictionary sliceTypeDict = d.subDict(transformKeys::SLICE_TYPE_DICT_KEY);
                singleTransformation.transformPoint = sliceTypeDict.lookup<point>(transformKeys::ORIGIN_KEY);
                singleTransformation.transformVector = sliceTypeDict.lookup<vector>(transformKeys::NORMAL_KEY);
                break; }
            case TransformationType::TRANSLATE:
                singleTransformation.transformCount = 0;
                singleTransformation.transformOffset = 0;
                singleTransformation.transformPoint = point(0,0,0);
                singleTransformation.transformVector = d.lookup<vector>(transformKeys::DISTANCE_KEY);
                break;
            case TransformationType::ROTATE:
                singleTransformation.transformCount = 0;
                singleTransformation.transformOffset = d.lookup<scalar>(transformKeys::ANGLE_DEGREES_KEY);
                singleTransformation.transformPoint = d.lookup<point>(transformKeys::ORIGIN_KEY);
                singleTransformation.transformVector = d.lookup<vector>(transformKeys::AXIS_KEY);
                break;
            case TransformationType::LINEAR_ARRAY:
                singleTransformation.transformCount = d.lookup<label>(transformKeys::COPY_COUNT_KEY);
                singleTransformation.transformOffset = d.lookup<scalar>(transformKeys::OFFSET_KEY);
                singleTransformation.transformPoint = point(0,0,0);
                singleTransformation.transformVector = d.lookup<vector>(transformKeys::DIRECTION_KEY);
                break;
            case TransformationType::SCALE:
                singleTransformation.transformCount = 0;
                singleTransformation.transformOffset = 0;
                singleTransformation.transformPoint = d.lookup<point>(transformKeys::ORIGIN_KEY);
                singleTransformation.transformVector = d.lookup<vector>(transformKeys::RATIO_KEY);
                break;
            }

            data_.transformations.append(singleTransformation);
        }
    }

    [[nodiscard]] const TransformObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new TransformDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] TransformInfo* copy() const override {
        return new TransformInfo(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        return isTypeAndDataEqualTo<TransformInfo>(other, data_);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
