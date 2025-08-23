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
    Foam::functionObjects::runTimeVis::StreamlineSourceInfo

Description
    Container for the streamline point source infos

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/objects/streamlineSources/streamlineSourceInfoFactory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class StreamlineSourceInfo
{
    std::shared_ptr<TemplatedStreamlineSourceInfo> ptr_;
public:
    StreamlineSourceInfo() = default;


    explicit StreamlineSourceInfo(const dictionary& streamlinesDict)
    {
        readDict(streamlinesDict);
    };

    void readDict(const dictionary& streamlinesDict)
    {
        ptr_.reset(StreamlineSourceInfoFactory::getFromDict(streamlinesDict));
    }

    void addToSources(std::vector<Id>& sources) const
    {
        ptr_->addToSources(sources);
    }

    [[nodiscard]] std::shared_ptr<StreamlineSourceProvider> createStreamlineSourceProvider() const
    {
        return ptr_->createStreamlineSourceProvider();
    }

    void computeAndAddToHash(size_t& hash) const
    {
        ptr_->computeAndAddToHash(hash);
    }

    bool operator==(const StreamlineSourceInfo& other) const
    {
        return ptr_->isDataEqualTo(other.ptr_.get());
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
