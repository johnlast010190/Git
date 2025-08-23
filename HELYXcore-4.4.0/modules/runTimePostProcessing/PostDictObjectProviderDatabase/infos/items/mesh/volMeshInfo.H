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
    Foam::functionObjects::runTimeVis::VolMeshInfo

Description
    Information for the volume mesh objects
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/mesh/meshInfo.H"
#include "itemDataSetProviders/mesh/volMeshDataSetProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class VolMeshInfo : public MeshInfo
{
public:
    VolMeshInfo
    (
        const string& name,
        const string& region,
        const dictionary& dict,
        const MeshAndFields& meshAndFields
    )
    :
        MeshInfo(name, region, dict, meshAndFields, Visualisation::solidColorDefaults(), ItemType::Value::PATCH)
    {};

    ~VolMeshInfo() override = default;

    [[nodiscard]] bool isInternalBoundary() const override {return false;};
    [[nodiscard]] bool isExternalBoundary() const override {return false;};
    [[nodiscard]] bool isProcessBoundary() const override {return false;};
    [[nodiscard]] bool isMeshBoundary() const override {return false;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new VolMeshDataSetProvider(id_.name, id_.region, meshAndFields_);
    };

    [[nodiscard]] VolMeshInfo* copy() const override {
        return new VolMeshInfo(*this);
    }

protected:
    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        auto casted = dynamic_cast<const VolMeshInfo*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return id_.name == casted->id_.name && id_.region == casted->id_.region;
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return true;};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
