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
    (c) 2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::MeridionalGridObjectData

Description
    Contains the dict data for the turbo background mesh/meridional plot
    post-processing object

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"
#include "meshes/primitiveShapes/point/point.H"

#include "baseClasses/id.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct MeridionalGridObjectData
{
    List<word> inletPatches;
    List<word> outletPatches;
    List<word> hubPatches;
    List<word> shroudPatches;

    label streamElements;
    label spanElements;

    point origin;
    vector axis;
    vector xDirection;

    void readDict(const dictionary &fieldSamplingDict)
    {
        inletPatches = fieldSamplingDict.lookup<List<word>>(turbopropKeys::INLET_PATCHES_KEY);
        outletPatches = fieldSamplingDict.lookup<List<word>>(turbopropKeys::OUTLET_PATCHES_KEY);
        hubPatches = fieldSamplingDict.lookup<List<word>>(turbopropKeys::HUB_PATCHES_KEY);
        shroudPatches = fieldSamplingDict.lookup<List<word>>(turbopropKeys::SHROUD_PATCHES_KEY);

        streamElements = fieldSamplingDict.lookup<label>(turbopropKeys::STREAM_ELEMENTS_KEY);
        spanElements = fieldSamplingDict.lookup<label>(turbopropKeys::SPAN_ELEMENTS_KEY);

        origin = fieldSamplingDict.lookup<point>(turbopropKeys::ORIGIN_KEY);
        axis = fieldSamplingDict.lookup<vector>(turbopropKeys::AXIS_KEY);
        xDirection = fieldSamplingDict.lookup<vector>(turbopropKeys::X_DIRECTION_KEY);
    }

    bool operator==(const MeridionalGridObjectData &other) const
    {
        // The patch sources are not equaled because they might refer to the same object even if different.
        // The objects pointed by the id should be compared instead.
        return streamElements == other.streamElements &&
               spanElements == other.spanElements &&
               origin == other.origin &&
               axis == other.axis &&
               xDirection == other.xDirection;
    }

    void computeAndAddToHash(size_t &hash) const
    {
        hasher::hash_combine(hash, streamElements);
        hasher::hash_combine(hash, spanElements);
        hasher::hash_combine(hash, origin);
        hasher::hash_combine(hash, axis);
        hasher::hash_combine(hash, xDirection);
        // The patch sources are hashed when we hash the sources
    }
};

} // End namespace Foam

// ************************************************************************* //
