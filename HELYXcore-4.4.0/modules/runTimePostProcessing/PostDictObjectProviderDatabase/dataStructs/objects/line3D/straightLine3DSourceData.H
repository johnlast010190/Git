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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::StraightLineSourceData

Description
    Contains the dict data necessary to define a straight line source for the line 3D

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "types/line3DSourceType.H"
#include "hash/hasher.H"

#include "primitives/Vector/vector/vector.H"
#include "db/dictionary/dictionary.H"

#include "engysLine3D.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct StraightLine3DSourceData
{
    label samplingPoints = 0;
    Vector<scalar> point1;
    Vector<scalar> point2;

    explicit StraightLine3DSourceData(const dictionary& line3DDict)
    {
        readDict(line3DDict);
    }

    void readDict(const dictionary& line3DDict)
    {
        samplingPoints = line3DDict.lookup<label>(line3DObjectKeys::N_SAMPLING_POINTS_KEY);
        point1 = line3DDict.lookup<Vector<scalar>>(line3DObjectKeys::POINT_1_KEY);
        point2 = line3DDict.lookup<Vector<scalar>>(line3DObjectKeys::POINT_2_KEY);
    }

    bool operator==(const StraightLine3DSourceData& other) const
    {
        return samplingPoints == other.samplingPoints &&
                       point1 == other.point1 &&
                       point2 == other.point2;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, Line3DSourceType::STRAIGHT_LINE);
        hasher::hash_combine(hash, samplingPoints);
        hasher::hash_combine(hash, point1);
        hasher::hash_combine(hash, point2);
    }
};

} // End namespace Foam

// ************************************************************************* //
