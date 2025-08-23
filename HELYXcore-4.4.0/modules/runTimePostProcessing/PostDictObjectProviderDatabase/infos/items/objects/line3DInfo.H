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
    Foam::functionObjects::runTimeVis::Line3DInfo

Description
    Information about the line 3D objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "infos/items/objects/line3DSources/line3DSourceInfoFactory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class Line3DInfo : public ObjectInfo
{
protected:

    std::shared_ptr<Line3DSourceInfo> sourceInfo_;
    const string& caseFolder_;

public:
    Line3DInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict,
        const string& caseFolder
    )
    :
        ObjectInfo(dict),
        caseFolder_(caseFolder)
    {
        Line3DInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~Line3DInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::LINE_3D);
        sourceInfo_.reset(Line3DSourceInfoFactory::getFromDict(dict, caseFolder_));
    }

    [[nodiscard]] const Line3DSourceInfo& getData() const { return *sourceInfo_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return sourceInfo_->createDataSetProvider(id_.name);
    };

    [[nodiscard]] Line3DInfo* copy() const override {
        return new Line3DInfo(*this);
    }

    [[nodiscard]] std::vector<Id> getSources() const override
    {
        std::vector<Id> sources;
        sources.push_back(source_);
        sourceInfo_->addExtraSources(sources);
        return sources;
    };

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        sourceInfo_->computeAndAddToHash(hash);
        return hash;
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        return isTypeAndDataEqualTo<Line3DInfo>(other, *sourceInfo_);
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        sourceInfo_->addItemRequirements(itemRequirements);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
