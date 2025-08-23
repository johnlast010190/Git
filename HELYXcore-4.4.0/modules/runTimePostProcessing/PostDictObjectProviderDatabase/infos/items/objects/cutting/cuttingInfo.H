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
    Foam::functionObjects::runTimeVis::CuttingInfo

Description
    Container for the cutting surface infos

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/cutting/cuttingInfoFactory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class CuttingInfo
{
    std::shared_ptr<TemplatedCuttingInfo> ptr_;
public:
    CuttingInfo() = default;


    CuttingInfo(const dictionary& sliceTypeDict, const string& caseFolder)
    {
        readDict(sliceTypeDict, caseFolder);
    };

    void readDict(const dictionary& sliceTypeDict, const string& caseFolder)
    {
        ptr_.reset(CuttingInfoFactory::getFromDict(sliceTypeDict, caseFolder));
    }

    [[nodiscard]] std::shared_ptr<CuttingSurfaceProvider> createCuttingSurfaceProvider() const
    {
        return ptr_->createCuttingSurfaceProvider();
    }

    void addItemRequirements(ItemRequirements& itemRequirements) const
    {
        ptr_->addItemRequirements(itemRequirements);
    }

    void computeAndAddToHash(size_t& hash) const
    {
        ptr_->computeAndAddToHash(hash);
    }

    bool operator==(const CuttingInfo& other) const
    {
        return ptr_->isDataEqualTo(other.ptr_.get());
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
