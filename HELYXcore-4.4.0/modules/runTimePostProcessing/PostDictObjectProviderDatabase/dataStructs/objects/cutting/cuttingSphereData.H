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
    Foam::functionObjects::runTimeVis::CuttingSphereData

Description
    Contains the dict data necessary to define a sphere slice

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "meshes/primitiveShapes/point/point.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"
#include "postDict/postDictKeys.H"
#include "types/cuttingType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct CuttingSphereData
{
    Vector<scalar> centre;
    scalar radius;

    explicit CuttingSphereData
    (
        const dictionary& sliceTypeDict
    )
    {
        centre = sliceTypeDict.lookup<point>(cuttingTypeKeys::CENTER_KEY);
        radius = sliceTypeDict.lookup<scalar>(cuttingTypeKeys::RADIUS_KEY);
    }

    bool operator==(const CuttingSphereData& other) const
    {
        return centre == other.centre &&
               radius == other.radius;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, CuttingType::SPHERE);
        hasher::hash_combine(hash, centre);
        hasher::hash_combine(hash, radius);
    }
};

} // End namespace

// ************************************************************************* //
