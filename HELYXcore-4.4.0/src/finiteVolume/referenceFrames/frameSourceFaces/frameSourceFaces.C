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

#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"
#include "containers/Lists/CompactListList/CompactListList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(frameSourceFaces, 0);
}


// * * * * * * * * * * * *  * * Constructors * * * * *  * * * * * * * * * * //

Foam::frameSourceFaces::frameSourceFaces
(
    const word& name,
    const fvMesh& mesh,
    const labelList& cells,
    const bool& active,
    const dictionary& dict
)
:
    regIOobject(IOobject(name, mesh.time().timeName(), mesh)),
    mesh_(mesh),
    cells_(cells),
    coeffsDict_(dict),
    nZoneFaces_(0),
    internalFaces_(mesh_.nFaces()),
    includedPatchFaces_(mesh_.boundaryMesh().size()),
    defaultMovingPatchFaces_(mesh_.boundaryMesh().size()),
    faceType_(mesh.nFaces(), notMoving),
    cellInFrame_(mesh_.nCells(), false)
{
    if (active)
    {
        UIndirectList<bool>(cellInFrame_, cells_) = true;
    }
    categorizeFaces();
    createFrameFaceLists();
}


void Foam::frameSourceFaces::categorizeFaces()
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (cellInFrame_[own[facei]] || cellInFrame_[nei[facei]])
        {
            faceType_[facei] = moving;
            nZoneFaces_++;
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!isA<emptyPolyPatch>(pp))
        {
            if
            (
                pp.coupled()
             || isA<cyclicAMIPolyPatch>(pp)
             || (
                   !includedPatchSet_.empty()
                && !includedPatchSet_.found(patchi)
                )
            )
            {
                forAll(pp, patchFacei)
                {
                    const label facei = pp.start() + patchFacei;
                    if (cellInFrame_[own[facei]])
                    {
                        faceType_[facei] = defaultMoving;
                        nZoneFaces_++;
                    }
                }
            }
            else
            {
                forAll(pp, patchFacei)
                {
                    const label facei = pp.start() + patchFacei;
                    if (cellInFrame_[own[facei]])
                    {
                        faceType_[facei] = moving;
                        nZoneFaces_++;
                    }
                }
            }
        }
    }
    syncTools::syncFaceList(mesh_, faceType_, maxEqOp<label>());
}


void Foam::frameSourceFaces::createFrameFaceLists()
{
    // Set internal faces
    label nInternal = 0;
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceType_[facei] == moving)
        {
            internalFaces_[nInternal++] = facei;
        }
    }
    internalFaces_.setSize(nInternal);

    // Construct lists of patch finite-volume faces in the frame
    {
        const fvBoundaryMesh& fvbMesh = mesh_.boundary();

        includedPatchFaces_.setSize(fvbMesh.size());
        defaultMovingPatchFaces_.setSize(fvbMesh.size());

        forAll(fvbMesh, patchi)
        {
            includedPatchFaces_[patchi].setSize(fvbMesh[patchi].size());
            defaultMovingPatchFaces_[patchi].setSize(fvbMesh[patchi].size());
        }

        labelList patchNFaces(fvbMesh.size(), 0);
        labelList defMovPatchNFaces(fvbMesh.size(), 0);

        for
        (
            label facei = mesh_.nInternalFaces(), bFacei = 0;
            facei < mesh_.nFaces();
            ++facei, ++bFacei
        )
        {
            if (faceType_[facei] == notMoving) continue;

            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            if (faceType_[facei] == moving)
            {
                forAll(patches, i)
                {
                    if (!fvbMesh[patches[i]].coupled())
                    {
                        includedPatchFaces_[patches[i]]
                            [patchNFaces[patches[i]]++] = patchFaces[i];
                    }
                    else
                    {
                        defaultMovingPatchFaces_[patches[i]]
                            [defMovPatchNFaces[patches[i]]++] = patchFaces[i];
                    }
                }
            }
            else if (faceType_[facei] == defaultMoving)
            {
                forAll(patches, i)
                {
                    defaultMovingPatchFaces_[patches[i]]
                        [defMovPatchNFaces[patches[i]]++] = patchFaces[i];
                }
            }
        }

        forAll(fvbMesh, patchi)
        {
            includedPatchFaces_[patchi].setSize(patchNFaces[patchi]);
            defaultMovingPatchFaces_[patchi].setSize(defMovPatchNFaces[patchi]);
        }
    }

    if (debug)
    {
        printDebugInformation();
    }
}


void Foam::frameSourceFaces::updateSourceFaces(const labelList& patchIDs)
{
    if (!patchIDs.empty())
    {
        includedPatchSet_ = patchIDs;
        categorizeFaces();
        createFrameFaceLists();
    }
}


const Foam::labelListList& Foam::frameSourceFaces::includedFaces() const
{
    return includedPatchFaces_;
}


const Foam::labelListList& Foam::frameSourceFaces::excludedFaces() const
{
    return defaultMovingPatchFaces_;
}


const Foam::labelList& Foam::frameSourceFaces::internalFaces() const
{
    return internalFaces_;
}


bool Foam::frameSourceFaces::writeData(Ostream& os) const
{
    return true;
}


void Foam::frameSourceFaces::printDebugInformation()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    faceSet internalFaces
    (
        mesh_,
        "internalFaces",
        labelHashSet(internalFaces_)
    );
    Pout<< "Writing " << internalFaces.size()
        << " internal faces in MRF zone to faceSet "
        << internalFaces.name() << endl;
    internalFaces.write();

    faceSet MRFFaces(mesh_, "includedFaces", 100);
    forAll(includedPatchFaces_, patchi)
    {
        forAll(includedPatchFaces_[patchi], i)
        {
            label patchFacei = includedPatchFaces_[patchi][i];
            MRFFaces.insert(patches[patchi].start() + patchFacei);
        }
    }
    Pout<< "Writing " << MRFFaces.size()
        << " patch faces in MRF zone to faceSet "
        << MRFFaces.name() << endl;
    MRFFaces.write();

    faceSet excludedFaces(mesh_, "excludedFaces", 100);
    forAll(defaultMovingPatchFaces_, patchi)
    {
        forAll(defaultMovingPatchFaces_[patchi], i)
        {
            label patchFacei = defaultMovingPatchFaces_[patchi][i];
            excludedFaces.insert(patches[patchi].start() + patchFacei);
        }
    }
    Pout<< "Writing " << excludedFaces.size()
        << " faces in MRF zone with special handling to faceSet "
        << excludedFaces.name() << endl;
    excludedFaces.write();
}


// ************************************************************************* //
