/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceSplitterObjectData

Description
    Contains the dict data for the testing of the surface angle splitter

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "baseClasses/id.H"
#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"
#include "hash/hasher.H"

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct SurfaceSplitterObjectData
{
    double angleDegrees;

    void readDict(const dictionary &dict)
    {
        angleDegrees = dict.lookup<double>(surfaceSplitterKeys::ANGLE_KEY);
    }

    bool operator==(const SurfaceSplitterObjectData &other) const
    {
        return angleDegrees == other.angleDegrees;
    }

    void computeAndAddToHash(size_t &hash) const
    {
        hasher::hash_combine(hash, angleDegrees);
    }
};

} // End namespace

// ************************************************************************* //
