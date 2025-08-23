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
    Foam::functionObjects::runTimeVis::SurfaceRotateData

Description
    Contains the dict data necessary to scale a surface

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#include "surfaceRotateData.H"
#include "postDict/postDictKeys.H"
#include "hash/hasher.H"
#include "vtkTransform.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{


void SurfaceRotateData::readDict(const dictionary &transformationDict)
{
    this->origin = transformationDict.lookup<Vector<scalar>>(surfaceTransformationKeys::ORIGIN_KEY);
    this->rotations = transformationDict.lookup<Vector<scalar>>(surfaceTransformationKeys::ROTATIONS_KEY);
}

bool SurfaceRotateData::operator==(const SurfaceTransformationData &other) const
{
    const auto* otherPtr = dynamic_cast<const SurfaceRotateData*>(&other);
    if (!otherPtr) return false;

    return this->origin == otherPtr->origin &&
           this->rotations == otherPtr->rotations;
}

void SurfaceRotateData::computeAndAddToHash(size_t& hash) const
{
    hasher::hash_combine(hash, origin);
    hasher::hash_combine(hash, rotations);
}

std::unique_ptr<SurfaceTransformationData> SurfaceRotateData::copy() const
{
    std::unique_ptr<SurfaceRotateData> unique = std::unique_ptr<SurfaceRotateData>(new SurfaceRotateData());
    unique->origin = origin;
    unique->rotations = rotations;
    return unique;
}

void SurfaceRotateData::addToVtkTransform(vtkTransform* transform) const
{
    transform->Translate(-origin.x(), -origin.y(), -origin.z());
    transform->RotateX(rotations.x());
    transform->RotateY(rotations.y());
    transform->RotateZ(rotations.z());
    transform->Translate(origin.x(), origin.y(), origin.z());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
