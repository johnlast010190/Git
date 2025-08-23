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
    Foam::functionObjects::runTimeVis::SurfaceCylinderData

Description
    Contains the dict data necessary for the cylinder surfaces

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "primitives/Vector/vector/vector.H"
#include "meshes/boundBox/boundBox.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

struct SurfacePlaneData
{
    Vector<scalar> basePoint;
    Vector<scalar> normalVector;
    boundBox meshBounds;
    scalar diagonal = NAN;

    void readDict(const dictionary &geometryDict, boundBox meshBounds_)
    {
        const dictionary &pointAndNormalDict = geometryDict.subDict(surfaceKeys::POINT_AND_NORMAL_DICT_KEY);
        basePoint = pointAndNormalDict.lookup<Vector<scalar>>(surfaceKeys::BASE_POINT_KEY);
        normalVector = pointAndNormalDict.lookup<Vector<scalar>>(surfaceKeys::NORMAL_VECTOR_KEY);
        diagonal = geometryDict.lookupOrDefault<scalar>(surfaceKeys::DIAGONAL_KEY, NAN);
        this->meshBounds = meshBounds_;
    }

    bool operator==(const SurfacePlaneData &other) const
    {
        return basePoint == other.basePoint &&
               normalVector == other.normalVector &&
               meshBounds == other.meshBounds &&
               std::isnan(diagonal) ? std::isnan(other.diagonal) : diagonal == other.diagonal;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, basePoint);
        hasher::hash_combine(hash, normalVector);
        if (!std::isnan(diagonal))
        {
            hasher::hash_combine(hash, diagonal);
        }
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
