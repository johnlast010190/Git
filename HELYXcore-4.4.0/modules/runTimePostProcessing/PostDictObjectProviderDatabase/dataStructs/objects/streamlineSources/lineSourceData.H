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
    Contains the dict data necessary to define a line source for the streamlines

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "types/streamlineSourceType.H"
#include "hash/hasher.H"

#include "primitives/Vector/vector/vector.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct LineSourceData
{
    label numberOfPoints;
    Vector<scalar> point1;
    Vector<scalar> point2;

    explicit LineSourceData(const dictionary& streamlinesDict)
    {
        numberOfPoints = streamlinesDict.lookup<label>(streamlinesSourceKeys::NUMBER_OF_POINTS_KEY);
        const dictionary& sourceDict = streamlinesDict.subDict(streamlinesSourceKeys::SOURCE_DICT_KEY);
        point1 = sourceDict.lookup<Vector<scalar>>(streamlinesSourceKeys::POINT_1_KEY);
        point2 = sourceDict.lookup<Vector<scalar>>(streamlinesSourceKeys::POINT_2_KEY);
    }

    bool operator==(const LineSourceData& other) const
    {
        return numberOfPoints == other.numberOfPoints &&
                       point1 == other.point1 &&
                       point2 == other.point2;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, StreamlineSourceType::LINE_SOURCE);
        hasher::hash_combine(hash, numberOfPoints);
        hasher::hash_combine(hash, point1);
        hasher::hash_combine(hash, point2);
    }
};

} // End namespace

// ************************************************************************* //
