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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "meshes/pointMesh/pointPatches/facePointPatch/facePointPatch.H"
#include "meshes/pointMesh/pointMesh.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"
#include "matrices/lduMatrix/lduAddressing/lduSchedule/lduSchedule.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointBoundaryMesh::pointBoundaryMesh
(
    const pointMesh& m,
    const polyBoundaryMesh& pbm
)
:
    pointPatchList(pbm.size()),
    mesh_(m)
{
    // Set boundary patches
    pointPatchList& Patches = *this;

    forAll(Patches, patchi)
    {
        Patches.set
        (
            patchi,
            facePointPatch::New(pbm[patchi], *this).ptr()
        );
    }
    // reset(pbm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointBoundaryMesh::findPatchID(const word& patchName) const
{
    return mesh()().boundaryMesh().findPatchID(patchName);
}


Foam::labelList Foam::pointBoundaryMesh::findIndices
(
    const keyType& key,
    const bool usePatchGroups
) const
{
    return mesh()().boundaryMesh().findIndices(key, usePatchGroups);
}


void Foam::pointBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initCalcGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initCalcGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }
}


void Foam::pointBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::pointBoundaryMesh::topoChange()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initTopoChange(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).topoChange(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initTopoChange(pBufs);
            }
            else
            {
                operator[](patchi).topoChange(pBufs);
            }
        }
    }
}


void Foam::pointBoundaryMesh::reset()
{
    const polyBoundaryMesh& boundaryMesh = mesh()().boundaryMesh();
    pointPatchList& Patches = *this;

    // Reset the number of patches in case the decomposition changed
    Patches.setSize(boundaryMesh.size());

    forAll(boundaryMesh, patchi)
    {
        // Construct new processor patches in case the decomposition changed
        if (isA<processorPolyPatch>(boundaryMesh[patchi]))
        {
            Patches.set
            (
                patchi,
                facePointPatch::New(boundaryMesh[patchi], *this).ptr()
            );
        }
    }
}


void Foam::pointBoundaryMesh::shuffle
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    pointPatchList::shuffle(newToOld);
    if (validBoundary)
    {
        topoChange();
    }
}


// ************************************************************************* //
