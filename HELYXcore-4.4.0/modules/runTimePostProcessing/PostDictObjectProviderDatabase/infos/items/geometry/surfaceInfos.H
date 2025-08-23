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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::TemplateSurfaceInfo

Description
    Information about the surface geometry objects, based on a template shared by them.
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"

#include "itemDataSetProviders/geometry/surfaceBoxDataSetProvider.H"
#include "itemDataSetProviders/geometry/surfaceCylinderDataSetProvider.H"
#include "itemDataSetProviders/geometry/surfaceFromFileDataSetProvider.H"
#include "itemDataSetProviders/geometry/surfacePlaneDataSetProvider.H"
#include "itemDataSetProviders/geometry/surfaceRingDataSetProvider.H"
#include "itemDataSetProviders/geometry/surfaceSphereDataSetProvider.H"

#include "dataStructs/geometry/surfaceTransformationsData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

template <typename Data, typename Provider>
class TemplateSurfaceInfo : public ItemInfo
{
protected:
    Data data_;
    SurfaceTransformationsData transformationsData_;

    TemplateSurfaceInfo() : ItemInfo(Visualisation::indexDefaults()) {};

    explicit TemplateSurfaceInfo
    (
        const dictionary& itemDict,
        const dictionary& geometryDict
    )
    :
        ItemInfo(itemDict, ItemType::Value::SURFACE, Visualisation::indexDefaults()),
        transformationsData_(geometryDict)
    {};

public:
    ~TemplateSurfaceInfo() override = default;

    const Data* getData() { return &data_;};

    [[nodiscard]] ItemDataSetProvider* createDataSetProvider() const override
    {
        return new Provider(id_.name, data_, transformationsData_);
    };

    [[nodiscard]] TemplateSurfaceInfo<Data, Provider>* copy() const override {
        return new TemplateSurfaceInfo<Data, Provider>(*this);
    }

protected:
    [[nodiscard]] size_t computeAndReturnSubclassItemInfoHash() const override
    {
        return data_.computeAndReturnHash();
    }

    bool isSubclassDataEqualTo(const ItemInfo* other) const override
    {
        auto casted = dynamic_cast<const TemplateSurfaceInfo<Data, Provider>*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };

    [[nodiscard]] bool needsGhostCellsWhenTransparent() const override {return false;};

    [[nodiscard]] bool needsGhostCellsWhenProfile() const override {return false;};
};

#define declareSurfaceInfoClass(name) class Surface##name##Info : public TemplateSurfaceInfo<Surface##name##Data, Surface##name##DataSetProvider>

declareSurfaceInfoClass(Box)
{
    public:
    SurfaceBoxInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(geometryDict);
    }
};

declareSurfaceInfoClass(Cylinder)
{
    public:
    SurfaceCylinderInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(geometryDict);
    }
};

declareSurfaceInfoClass(FromFile)
{
    public:
    SurfaceFromFileInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict,
        const string& caseFolder
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(itemDict, geometryDict, caseFolder);
    }
};

declareSurfaceInfoClass(Plane)
{
    public:
    SurfacePlaneInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict,
        const boundBox& meshBounds
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(geometryDict, meshBounds);
    }
};

declareSurfaceInfoClass(Sphere)
{
    public:
    SurfaceSphereInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(geometryDict);
    }
};

declareSurfaceInfoClass(Ring)
{
    public:
    SurfaceRingInfo
    (

        const dictionary& itemDict,
        const dictionary& geometryDict
    )
    :
        TemplateSurfaceInfo(itemDict, geometryDict)
    {
        data_.readDict(geometryDict);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
