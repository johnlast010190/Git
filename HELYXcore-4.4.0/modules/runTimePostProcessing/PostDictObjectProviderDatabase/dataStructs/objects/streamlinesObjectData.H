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
    Foam::functionObjects::runTimeVis::StreamlinesObjectData

Description
    Contains the dict data for the streamlines, expect for the source, which is
     contained in one of the streamline source objects

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct StreamlinesObjectData
{
    foamField vectorField {};
    scalar maxLength {};
    label maxSteps {};
    scalar radius {};

    void readDict(const dictionary& streamlinesDict)
    {
        vectorField = foamField(streamlinesDict.lookup<word>(streamlinesObjectKeys::VECTOR_FIELD_KEY));
        maxLength = streamlinesDict.lookup<scalar>(streamlinesObjectKeys::MAX_LENGTH_KEY);
        maxSteps = streamlinesDict.lookup<label>(streamlinesObjectKeys::MAX_STEPS_KEY);
        radius = streamlinesDict.lookup<scalar>(streamlinesObjectKeys::RADIUS_KEY);
    }

    bool operator==(const StreamlinesObjectData& other) const
    {
        return vectorField == other.vectorField &&
                 maxLength == other.maxLength &&
                  maxSteps == other.maxSteps &&
                    radius == other.radius;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, vectorField);
        hasher::hash_combine(hash, maxLength);
        hasher::hash_combine(hash, maxSteps);
        hasher::hash_combine(hash, radius);
    }
};

} // End namespace

// ************************************************************************* //
