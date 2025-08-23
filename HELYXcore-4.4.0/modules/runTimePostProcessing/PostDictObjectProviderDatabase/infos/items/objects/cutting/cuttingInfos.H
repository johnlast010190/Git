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
    Foam::functionObjects::runTimeVis::TemplatedCuttingInfo

Description
    Contains the cutting infos, based on a template class

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/objects/cutting/cuttingSurfaceProviders.H"
#include "baseClasses/itemRequirements.H"
#include "db/dictionary/dictionary.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class TemplatedCuttingInfo
{
public:
    virtual ~TemplatedCuttingInfo() = default;


    // Member Functions

    [[nodiscard]] virtual std::shared_ptr<CuttingSurfaceProvider> createCuttingSurfaceProvider() const = 0;

    virtual void computeAndAddToHash(size_t& hash) const  = 0;

    virtual bool isDataEqualTo(const TemplatedCuttingInfo* other) const = 0;

    virtual void addItemRequirements(ItemRequirements& itemRequirements) const {};
};

template <typename Data, typename Provider, typename ...Arg>
class TemplateCuttingInfo : public TemplatedCuttingInfo
{
private:
    Data data_;

public:
    __attribute__((unused)) explicit TemplateCuttingInfo(const dictionary& cuttingTypeDict, const Arg& ...arg)
    :
        data_(cuttingTypeDict, arg...)
    {};

    [[nodiscard]] std::shared_ptr<CuttingSurfaceProvider> createCuttingSurfaceProvider() const override
    {
        return std::shared_ptr<CuttingSurfaceProvider>(new Provider(data_));
    };

    void computeAndAddToHash(size_t& hash) const override
    {
        data_.computeAndAddToHash(hash);
    }

    bool isDataEqualTo(const TemplatedCuttingInfo* other) const override
    {
        auto casted = dynamic_cast<const TemplateCuttingInfo<Data, Provider, Arg...>*>(other);
        if (casted == nullptr)
        {
            return false;
        }
        return data_ == casted->data_;
    };
};

// Declaration of the individual cutting infos
typedef TemplateCuttingInfo<CuttingPlaneData, CuttingPlaneProvider> CuttingPlaneInfo;
typedef TemplateCuttingInfo<CuttingBoxData, CuttingBoxProvider> CuttingBoxInfo;
typedef TemplateCuttingInfo<CuttingCylinderData, CuttingCylinderProvider> CuttingCylinderInfo;
typedef TemplateCuttingInfo<CuttingSphereData, CuttingSphereProvider> CuttingSphereInfo;

class CuttingFileInfo : public TemplateCuttingInfo<CuttingFileData, CuttingFileProvider, string>
{
public:
    __attribute__((unused)) explicit CuttingFileInfo(const dictionary& cuttingTypeDict, const string& arg)
        :
        TemplateCuttingInfo<CuttingFileData, CuttingFileProvider, string>(cuttingTypeDict, arg)
    {};

    void addItemRequirements(ItemRequirements& itemRequirements) const override
    {
        // needs it because the method actually depends on a probing of dataset point values, where the points
        // can and often do fall right on the process boundaries, which could be a miss if there were no
        // ghost cells
        itemRequirements.setNeedsGhostCellsTrue();
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
