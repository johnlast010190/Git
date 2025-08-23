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
    Foam::functionObjects::runTimeVis::GlyphInfo

Description
    Information about the glyph objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/glyphDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class GlyphInfo : public ObjectInfo
{
protected:

    GlyphInfo() = default;
    GlyphObjectData data_;

public:
    GlyphInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        GlyphInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~GlyphInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::GLYPH);
        data_.readDict(dict);
    }

    [[nodiscard]] const GlyphObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new GlyphDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] GlyphInfo* copy() const override {
        return new GlyphInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        return isTypeAndDataEqualTo<GlyphInfo>(other, data_);
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        foamField orientationField = foamField::withPointAssociation(data_.orientationField);
        itemRequirements.addToRequiredFields(orientationField);
        if (data_.enableScalingByField)
        {
            foamField scaleField = foamField::withPointAssociation(data_.scalingField);
            itemRequirements.addToRequiredFields(scaleField);
        }
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
