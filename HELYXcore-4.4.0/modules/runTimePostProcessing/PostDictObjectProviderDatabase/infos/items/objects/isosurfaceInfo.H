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
    Foam::functionObjects::runTimeVis::IsoSurfaceInfo

Description
    Information about the iso surface objects.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/objectInfo.H"
#include "itemDataSetProviders/objects/isosurfaceDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class IsoSurfaceInfo : public ObjectInfo
{
protected:

    IsoSurfaceInfo() = default;
    IsoSurfaceObjectData data_;

public:
    IsoSurfaceInfo
    (
        const dictionary& dict,
        const dictionary& objectVisDict
    )
    :
        ObjectInfo(dict)
    {
        IsoSurfaceInfo::updateObjectData(dict);
        readVisualisationDict(objectVisDict);
    };
    ~IsoSurfaceInfo() override = default;

    void updateObjectData(const dictionary& dict) override
    {
        checkObjectType(ItemType::Value::ISO_SURFACE);

        data_.values = dict.lookup<List<scalar>>("values");
        data_.scalarField = foamField(dict.lookup<word>("scalarField"));
        data_.enableSmoothing = dict.lookup<bool>("enableSmoothing");
        data_.iterations = dict.lookup<label>("iterations");
        data_.convergence = dict.lookup<scalar>("convergence");
    }

    [[nodiscard]] const IsoSurfaceObjectData& getData() const { return data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new IsoSurfaceDataSetProvider(id_.name, data_);
    };

    [[nodiscard]] IsoSurfaceInfo* copy() const override {
        return new IsoSurfaceInfo(*this);
    }

protected:

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        return isTypeAndDataEqualTo<IsoSurfaceInfo>(other, data_);
    };

    void addExtraItemRequirements(ItemRequirements& itemRequirements) const override
    {
        foamField scalarField = foamField::withPointAssociation(data_.scalarField);
        itemRequirements.addToRequiredFields(scalarField);
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
