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
    Foam::functionObjects::runTimeVis::LineSourceData

Description
    Contains the dict data necessary to define a patch source for the streamlines

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "types/streamlineSourceType.H"
#include "baseClasses/id.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct PatchSourceData
{
    label numberOfPoints;
    Id patchId;

    explicit PatchSourceData(const dictionary& streamlinesDict)
    {
        numberOfPoints = streamlinesDict.lookup<label>(streamlinesSourceKeys::NUMBER_OF_POINTS_KEY);
        const dictionary& sourceDict = streamlinesDict.subDict(streamlinesSourceKeys::SOURCE_DICT_KEY);
        const dictionary& patchIdDict = sourceDict.subDict(streamlinesSourceKeys::PATCH_ID_KEY);
        patchId = Id(patchIdDict);
    }

    bool operator==(const PatchSourceData& other) const
    {
        // The patch source is not equaled because they might refer to the same object even if different.
        // The objects pointed by the id should be compared instead.
        return numberOfPoints == other.numberOfPoints;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, StreamlineSourceType::PATCH_SOURCE);
        hasher::hash_combine(hash, numberOfPoints);
        // The patch source is hashed when we hash the sources
    }
};

} // End namespace

// ************************************************************************* //
