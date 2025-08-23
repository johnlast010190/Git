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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrame.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(coordinateFrame, 0);
    defineRunTimeSelectionTable(coordinateFrame, dictionary);
    addToRunTimeSelectionTable(coordinateFrame, coordinateFrame, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateFrame::coordinateFrame
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& frameName
)
:
    MeshObject<fvMesh, UpdateableMeshObject, coordinateFrame>(mesh, frameName),
    coordinateFrameState(mesh.thisDb(), dict, frameName),
    coorFrameReg_(mesh),
    frameDict_(dict),
    parentFrameName_
    (
        frameDict_.lookupOrAddDefault<word>("parentFrameName", word::null)
    )
{
    if (anyIncremental())
    {
        isIncrementalMotion() = true;
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::coordinateFrame::~coordinateFrame()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordinateFrame::initParents() const
{
    if (parents_.empty())
    {
        UPtrList<coordinateFrame> parents;
        label nNestedFrames = 0;
        bool nextFrame = true;
        coordinateFrame* framePtr = const_cast<coordinateFrame*>(this);
        while (nextFrame)
        {
            ++nNestedFrames;
            parents.setSize(nNestedFrames);
            parents.set(nNestedFrames - 1, const_cast<coordinateFrame*>(framePtr));
            if (framePtr->validParentFrame())
            {
                framePtr =
                    const_cast<coordinateFrame*>(&framePtr->parentFrame());
            }
            else
            {
                nextFrame = false;
            }
        }

        parents_.setSize(nNestedFrames);
        label n = 0;
        for (label i = parents.size() - 1; i >= 0; --i)
        {
            parents_.set(n, &parents[i]);
            ++n;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateFrame::localCoordSystem() const
{
    if
    (
        validParentFrame()
     && coorSys().coordinateTransforms().isLocalOnRead()
    && !coorSys().coordinateTransforms().isLocal()
    )
    {
        coordinateSystem localCS(parents().last().coorSys().clone());
        const label parentIndex = parents().size() - 2;

        if (parents()[parentIndex].coorSys().coordinateTransforms().isLocal())
        {
            FatalErrorInFunction
                << "Assumes that all frames are in global coordinates"
                << nl << abort(FatalError);
        }

        frameData::transformToLocal(localCS, parents()[parentIndex].coorSys());

        return localCS.clone();
    }
    else
    {
        return coorSys().clone();
    }
}


bool Foam::coordinateFrame::validParentFrame() const
{
    if (mesh_.foundObject<coordinateFrame>(parentFrameName_))
    {
        parentFramePtr_ =
            mesh_.lookupObjectRefPtr<coordinateFrame>(parentFrameName_);
        return true;
    }
    else if (parentFrameName_ != word::null)
    {
        parentFramePtr_ = &coordinateFrame::New(mesh_, parentFrameName_);
        return true;
    }
    return false;
}


const Foam::septernion& Foam::coordinateFrame::transformation
(
    bool addParentFrames,
    label nCorr
) const
{
    globalTransform_ = septernion::I;

    if (addParentFrames)
    {
        forAll(parents(), i)
        {
            if (parents()[i].isDynamic())
            {
                globalTransform_ *=
                    parents()[i].decoupledTransformation(nCorr);
            }
        }
    }
    else
    {
        if (isDynamic())
        {
            globalTransform_ = decoupledTransformation(nCorr);
        }
    }

    return globalTransform_;
}


void Foam::coordinateFrame::updateState(bool construction) const
{
    if (validParentFrame())
    {
        parentFrame().updateState(construction);
    }

    if (!isUpdated() && validParentFrame())
    {
        // Request construction of old times
        oldTime().oldTime();

        // First store old time data, then update
        storeOldTimes();

        // Frame can get updated based on parent frame
        if (!construction)
        {
            updateCoordinateSystem();
        }

        // Update transformation, velocities and accelerations
        const label nCorr = nFrameCorrector() + 1;
        transformations_ = List<septernion>(nCorr, septernion::I);
        velocities_ = List<vectorTuple>(nCorr, vectorTuple(Zero, Zero));
        accelerations_ = List<vectorTuple>(nCorr, vectorTuple(Zero, Zero));

        updateIndex_ = obr_.time().timeIndex();
        if (obr_.time().writeTime())
        {
            coordinateFrameState::write();
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::coordinateFrame::frameVelocity
(
    const volVectorField& positions,
    bool addParentFrames
) const
{
    updateState();

    // Frame velocity
    tmp<volVectorField> U
    (
        volVectorField::New("Urf", mesh_, dimensionedVector(dimVelocity, Zero))
    );

    if (addParentFrame(addParentFrames))
    {
        U.ref() += parentFrame().frameVelocity(positions, addParentFrames);
    }

    return U;
}


Foam::tmp<Foam::vectorField> Foam::coordinateFrame::frameVelocity
(
    const vectorField& positions,
    bool addParentFrames
) const
{
    updateState();

    tmp<vectorField> U(new vectorField(positions.size(), Zero));
    if (addParentFrame(addParentFrames))
    {
        U.ref() += parentFrame().frameVelocity(positions, addParentFrames);
    }

    return U;
}


Foam::vector Foam::coordinateFrame::frameVelocity
(
    const vector& position,
    bool addParentFrames
) const
{
    updateState();

    if (addParentFrame(addParentFrames))
    {
        return parentFrame().frameVelocity(position, addParentFrames);
    }

    return Zero;
}


void Foam::coordinateFrame::attachPatch(const label& patchi) const
{
    coorFrameReg_.attachPatch(patchi);

    if (validParentFrame())
    {
        parentFrame().attachPatch(patchi);
    }
}


void Foam::coordinateFrame::inletFluxVelocity
(
    vectorField& Up,
    const label patchi,
    const objectRegistry& obr,
    const word& UName,
    bool addedMeshPhi
) const
{
    const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];

    // For nested dynamic frames remove normal component of velocity for all
    // frames but add meshPhi normal velocity component only ones
    if (mesh_.moving() && isDynamic())
    {
        const vectorField& nf = mesh_.boundaryMesh()[patchi].faceNormals();
        tmp<vectorField> Uframe = frameVelocity(Cf, false);
        Uframe.ref() -= nf*(nf & Uframe());

        // Mesh moving - normal component of velocity is removed and replaced
        // with the mesh flux velocity
        if (!addedMeshPhi)
        {
            const volVectorField& U = obr.lookupObject<volVectorField>(UName);
            const fvPatch& pp = U.boundaryField()[patchi].patch();
            const scalarField phip(fvc::meshPhi(U, pp.index()));

            tmp<scalarField> Un = phip/(pp.magSf() + VSMALL);
            Uframe.ref() += nf*Un;
            Up += Uframe;
            addedMeshPhi = true;
        }
    }
    else
    {
        Up += frameVelocity(Cf, false);
    }

    if (validParentFrame())
    {
        parentFramePtr_->inletFluxVelocity
        (
            Up,
            patchi,
            obr,
            UName,
            addedMeshPhi
        );
    }
}


void Foam::coordinateFrame::calculateBoundaryVelocity
(
    vectorField& Up,
    const label patchi,
    const objectRegistry& obr,
    const word& UName,
    const bool inletFlux
) const
{
    const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];
    if (inletFlux)
    {
        inletFluxVelocity(Up, patchi, obr, UName);
        return;
    }
    else if (!isAttachToMRF(patchi))
    {
        attachPatch(patchi);
    }
    const polyPatch& pp = mesh_.boundaryMesh()[patchi];
    vectorField Uframe(frameVelocity(Cf, false));
    boolList isFaceInMRF(Up.size(), false);

    forAll(coorFrameReg_.registeredNames(), objI)
    {
        const word& objName = coorFrameReg_.registeredNames()[objI];
        if (mesh_.foundObject<frameSourceFaces>(objName))
        {
            const frameSourceFaces& frameSrcFaces =
                mesh_.lookupObject<frameSourceFaces>(objName);
            const labelList& patchFaces = frameSrcFaces.includedFaces()[patchi];
            forAll(patchFaces, i)
            {
                const label patchFacei = patchFaces[i];
                Up[patchFacei] += Uframe[patchFacei];
                isFaceInMRF[patchFacei] = true;
            }
        }
    }
    const vectorField& nfs = mesh_.boundaryMesh()[patchi].faceNormals();

    forAll(pp, patchFacei)
    {
        if (!isFaceInMRF[patchFacei])
        {
            const vector& nf = nfs[patchFacei];
            Up[patchFacei] +=
                Uframe[patchFacei] - nf*(nf & Uframe[patchFacei]);
        }
    }

    // If the parent frame isn't valid then it is the last frame in the chain
    // and the mesh flux velocity is added
    if (validParentFrame())
    {
        parentFrame().calculateBoundaryVelocity
        (
            Up,
            patchi,
            obr,
            UName,
            inletFlux
        );
    }
    else if (mesh_.moving())
    {
        // Mesh moving - normal component of velocity is removed and replaced
        // with the mesh flux velocity
        const volVectorField& U = obr.lookupObject<volVectorField>(UName);
        const fvPatch& pp = U.boundaryField()[patchi].patch();
        const scalarField phip(fvc::meshPhi(U, pp.index()));

        tmp<scalarField> Un = phip/(pp.magSf() + VSMALL);
        Up += nfs*Un;
    }
}


bool Foam::coordinateFrame::isAttachToMRF(const label& patchi) const
{
    if (coorFrameReg_.isAttachToMRF(patchi))
    {
        return true;
    }
    if (validParentFrame())
    {
        return parentFrame().isAttachToMRF(patchi);
    }
    return false;
}


// ************************************************************************* //
