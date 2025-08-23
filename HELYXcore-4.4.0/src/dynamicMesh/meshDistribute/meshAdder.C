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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "meshDistribute/meshAdder.H"
#include "polyMeshAdder/faceCoupleInfo.H"
#include "fvMesh/fvMesh.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(meshAdder, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::meshAdder::calcPatchMap
(
    const label oldStart,
    const label oldSize,
    const labelList& oldToNew,
    const polyPatch& newPatch,
    const label unmappedValue
)
{
    labelList newToOld(newPatch.size(), unmappedValue);

    label newStart = newPatch.start();
    label newSize = newPatch.size();

    for (label i = 0; i < oldSize; i++)
    {
        label newFacei = oldToNew[oldStart+i];

        if (newFacei >= newStart && newFacei < newStart+newSize)
        {
            newToOld[newFacei-newStart] = i;
        }
    }
    return newToOld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::meshAdder::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary,
    const bool fullyMapped
)
{
    // Store old mesh0 point maps
    labelListList oldMeshPoints0;
    const bool havePointMesh =
        mesh0.foundObject<pointMesh>(pointMesh::typeName);
    if (havePointMesh)
    {
        const polyBoundaryMesh& pbm0 = mesh0.boundaryMesh();
        oldMeshPoints0.setSize(pbm0.size());
        forAll(pbm0, patchi)
        {
            oldMeshPoints0[patchi] = pbm0[patchi].meshPoints();
        }
    }

    // Resulting merged mesh (polyMesh only!)
    autoPtr<mapAddedPolyMesh> mapPtr
    (
        polyMeshAdder::add
        (
            mesh0,
            mesh1,
            coupleInfo,
            validBoundary
        )
    );

    // Adjust the fvMesh part.
    const polyBoundaryMesh& patches = mesh0.boundaryMesh();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh0.boundary());
    fvPatches.setSize(patches.size());
    forAll(patches, patchi)
    {
        fvPatches.set(patchi, fvPatch::New(patches[patchi], fvPatches));
    }

    if (havePointMesh)
    {
        // Recreate point mesh
        const pointMesh& pointMesh0 = pointMesh::New(mesh0);

        MapPointFields<scalar>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<vector>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<sphericalTensor>
        (
            mapPtr,
            pointMesh0,
            oldMeshPoints0,
            mesh1
        );
        MapPointFields<symmTensor>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
        MapPointFields<tensor>(mapPtr, pointMesh0, oldMeshPoints0, mesh1);
    }

    // clear the geometry data so re-computation gets triggered during mapping
    mesh0.clearOut();

    // Do the mapping of the stored fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MapVolFields<scalar>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<vector>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<sphericalTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<symmTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<tensor>(mapPtr, mesh0, mesh1, fullyMapped);

    MapVolFields<vector1>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<vector4>(mapPtr, mesh0, mesh1, fullyMapped);
    MapVolFields<tensor4>(mapPtr, mesh0, mesh1, fullyMapped);


    MapSurfaceFields<scalar>(mapPtr, mesh0, mesh1, fullyMapped);
    MapSurfaceFields<vector>(mapPtr, mesh0, mesh1, fullyMapped);
    MapSurfaceFields<sphericalTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    MapSurfaceFields<symmTensor>(mapPtr, mesh0, mesh1, fullyMapped);
    MapSurfaceFields<tensor>(mapPtr, mesh0, mesh1, fullyMapped);


    MapDimFields<scalar>(mapPtr, mesh0, mesh1);
    MapDimFields<vector>(mapPtr, mesh0, mesh1);
    MapDimFields<sphericalTensor>(mapPtr, mesh0, mesh1);
    MapDimFields<symmTensor>(mapPtr, mesh0, mesh1);
    MapDimFields<tensor>(mapPtr, mesh0, mesh1);

    MapDimFields<vector1>(mapPtr, mesh0, mesh1);
    MapDimFields<vector4>(mapPtr, mesh0, mesh1);
    MapDimFields<tensor4>(mapPtr, mesh0, mesh1);

    return mapPtr;
}


// ************************************************************************* //
