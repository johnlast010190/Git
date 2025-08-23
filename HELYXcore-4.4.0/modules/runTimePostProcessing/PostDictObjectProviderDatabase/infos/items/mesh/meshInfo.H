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
    Foam::functionObjects::runTimeVis::MeshInfo

Description
    Base class for the mesh item infos
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "baseClasses/meshAndFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class MeshInfo : public ItemInfo
{
protected:

    MeshAndFields meshAndFields_;

    MeshInfo
    (
        const string& name,
        const string& region,
        const dictionary& dict,
        const MeshAndFields& meshAndFields,
        const Visualisation* defaultVisualisation,
        ItemType::Value type
    )
    :
        MeshInfo(name, region, meshAndFields, defaultVisualisation, type)
    {
        visualisation_.readFromDict(dict);
    };

    MeshInfo
    (
        const string& name,
        const string& region,
        const MeshAndFields& meshAndFields,
        const Visualisation* defaultVisualisation,
        ItemType::Value type
    )
    :
        ItemInfo(defaultVisualisation),
        meshAndFields_(meshAndFields)
    {
        id_.name = name;
        id_.region = region;
        id_.type = type;
    };

    ~MeshInfo() override = default;

    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        size_t hash = 0;
        hasher::hash_combine(hash, id_.name);
        hasher::hash_combine(hash, id_.region);
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace

// ************************************************************************* //
