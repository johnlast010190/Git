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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrameRegistry/coordinateFrameRegistry.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateFrameRegistry::coordinateFrameRegistry(const fvMesh& mesh)
:
    mesh_(mesh),
    registeredNames_(),
    registeredPatches_(),
    faceSetNeedUpdate_(false)
{}


void Foam::coordinateFrameRegistry::attachObject(const word& name) const
{
    registeredNames_.append(name);
    faceSetNeedUpdate_ = true;
    updateFaces();
}


bool Foam::coordinateFrameRegistry::attachPatch(label patchi) const
{
    bool newPatch = true;
    faceSetNeedUpdate_ = false;
    forAll(registeredPatches_, pI)
    {
        if (registeredPatches_[pI] == patchi)
        {
            newPatch = false;
        }
    }
    if (newPatch)
    {
        registeredPatches_.append(patchi);
        faceSetNeedUpdate_ = true;
    }
    updateFaces();

    return faceSetNeedUpdate_;
}


bool Foam::coordinateFrameRegistry::isAttachToMRF(label patchi) const
{
    forAll(registeredNames_, i)
    {
        const word& name = registeredNames_[i];
        if
        (
            mesh_.foundObject<frameSourceFaces>(name)
         && mesh_.lookupObject<frameSourceFaces>
            (
                name
            ).includedFaces()[patchi].size()
        )
        {
            return true;
        }
    }
    return false;
}


void Foam::coordinateFrameRegistry::updateFaces() const
{
    if (faceSetNeedUpdate_)
    {
        forAll(registeredNames_, i)
        {
            const word& name = registeredNames_[i];
            if (mesh_.foundObject<frameSourceFaces>(name))
            {
                mesh_.lookupObjectRef<frameSourceFaces>
                (
                    name
                ).updateSourceFaces(registeredPatches_);
            }
        }
    }
}


// ************************************************************************* //
