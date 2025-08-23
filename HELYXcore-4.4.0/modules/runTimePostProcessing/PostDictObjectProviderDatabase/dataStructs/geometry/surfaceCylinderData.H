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
    Foam::functionObjects::runTimeVis::SurfaceCylinderData

Description
    Contains the dict data necessary for the cylinder surfaces

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "primitives/Vector/vector/vector.H"
#include "hash/hasher.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

struct SurfaceCylinderData {
    Vector<scalar> point1;
    Vector<scalar> point2;
    scalar radius = 0;

    void readDict(const dictionary& geometryDict)
    {
        point1 = geometryDict.lookup<Vector<scalar>>(surfaceKeys::POINT_1_KEY);
        point2 = geometryDict.lookup<Vector<scalar>>(surfaceKeys::POINT_2_KEY);
        radius = geometryDict.lookup<scalar>(surfaceKeys::RADIUS_KEY);
    }

    bool operator==(const SurfaceCylinderData& other) const
    {
        return point1 == other.point1 &&
                point2 == other.point2 &&
                radius == other.radius;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, point1);
        hasher::hash_combine(hash, point2);
        hasher::hash_combine(hash, radius);
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
