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
    (c) 2014-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "layerManipulate/layerManipulate.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"
#include "algorithms/PointEdgeWave/PointEdgeWave.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "hessianMeshOptimization/leastSquaresCurveFit/leastSquaresCurveFit.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "edgeClassification/edgeClassification.H"
#include "meshes/polyMesh/polyMeshCheck/polyMeshTools.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

namespace Foam
{

defineTypeNameAndDebug(layerManipulate, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeFaceZonePatch()
{
    const fvMesh& mesh = meshRefiner_.mesh();
    labelList addressing;

    // Count faces.
    label nFaces = 0;
    forAll(mesh.faceZones(), zoneID)
    {
        nFaces += mesh.faceZones()[zoneID].size();
    }
    addressing.setSize(nFaces);

    nFaces = 0;
    forAll(mesh.faceZones(), zoneID)
    {
        const labelList& zoneFaces = mesh.faceZones()[zoneID];
        forAll(zoneFaces, i)
        {
            label facei = zoneFaces[i];
            addressing[nFaces++] = facei;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), addressing),
            mesh.points()
        )
    );
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeLayerPatch()
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    boolList internalBafflePatches(patches.size(), false);
    forAll(mesh.faceZones(), zonei)
    {
        label mpi, spi;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = meshRefiner_.getFaceZoneInfo
        (
            mesh.faceZones()[zonei].name(),
            mpi,
            spi,
            fzType
        );
        if
        (
            hasInfo &&
            (
                fzType == surfaceZonesInfo::INTERNAL
                || fzType == surfaceZonesInfo::BAFFLE
            )
        )
        {
            internalBafflePatches[mpi] = true;
            internalBafflePatches[spi] = true;
        }
    }

    const labelList& patchToNLayers = layerParams_.numLayers();

    DynamicList<label> adaptPatches(patches.size());

    forAll(patches, patchI)
    {
        if
        (
            !patches[patchI].coupled()
            && patchToNLayers[patchI] != -1
            && !internalBafflePatches[patchI]
        )
        {
            adaptPatches.append(patchI);
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatches
        )
    );
}

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeGrownUpPatch(const bool onlySnap)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& patchToNLayers = layerParams_.numLayers();
    const List<Switch>& reSnap = layerParams_.reSnap();

    DynamicList<label> grownUpPatches(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && patchToNLayers[patchI] == -1)
        {
            if (onlySnap)
            {
                if (reSnap[patchI])
                {
                    grownUpPatches.append(patchI);
                }
            }
            else
            {
                grownUpPatches.append(patchI);
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        meshRefinement::makePatch
        (
            mesh,
            grownUpPatches
        )
    );
}

// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::layerManipulate::isMaxEdge
(
    const fvMesh& mesh,
    const List<pointData>& pointWallDist,
    const label edgei,
    const scalar minCos
)
{
    int dummyTrackData = 0;
    const edge& e = mesh.edges()[edgei];
    if
    (
        !pointWallDist[e[0]].valid(dummyTrackData)
        || !pointWallDist[e[1]].valid(dummyTrackData)
    )
    {
        return false;
    }

    const pointField& points = mesh.points();

    // Do not mark edges with one side on moving wall.

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < VSMALL)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < VSMALL)
    {
        return false;
    }

    //- Detect based on extrusion vector differing for both endpoints
    //  the idea is that e.g. a sawtooth wall can still be extruded
    //  successfully as long as it is done all to the same direction.
    if ((pointWallDist[e[0]].v() & pointWallDist[e[1]].v()) < minCos)
    {
        vector origVec =
            pointWallDist[e[1]].origin() - pointWallDist[e[0]].origin();
        scalar magOrigVec = mag(origVec);

        if (magOrigVec > VSMALL)
        {
            origVec /= magOrigVec;
            if
            (
                (origVec & pointWallDist[e[1]].v()) > 0
                && (origVec & pointWallDist[e[0]].v()) < 0
             )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::layerManipulate::calculateLayerStacks
(
    const indirectPrimitivePatch& pp,
    const boolList& boundaryFaces,
    const boolList& boundaryEdges,
    labelList& cellLayerNumber,
    labelList& layerFaceType,
    labelList& layerPointType,
    labelList& stackEdgeLayer
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();

    const labelIOList& layerCells = meshRefiner_.layerCells()();

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] == -1)
        {
            cellLayerNumber[celli] = -1;
        }
    }

    boolList bPoints(mesh.nPoints(), false);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        if (layerCells[owners[meshFaceI]] < 0)
        {
            continue;
        }

        layerFaceType[meshFaceI] = 2;

        label celli = owners[meshFaceI];

        if (cellLayerNumber[celli] == 0)
        {
            cellLayerNumber[celli]++;
        }

        face f = pp[i];
        forAll(f, fp)
        {
            label meshPointI = f[fp];
            bPoints[meshPointI] = true;
            const labelList& pEdges = mesh.pointEdges()[meshPointI];

            forAll(pEdges, pEI)
            {
                labelHashSet cEdges(mesh.cellEdges()[celli]);
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    if
                    (
                        !boundaryEdges[edgei] && stackEdgeLayer[edgei] == 0
                        && cEdges.found(edgei)
                     )
                    {
                        stackEdgeLayer[edgei]++;
                    }
                }
            }
        }
        const labelList& fEdges = mesh.faceEdges()[meshFaceI];
        labelHashSet cFaces(mesh.cells()[celli]);
        forAll(fEdges, fEI)
        {
            label edgei = fEdges[fEI];

            const labelList& eFaces = mesh.edgeFaces()[edgei];
            forAll(eFaces, eFI)
            {
                label facei = eFaces[eFI];
                if
                (
                    facei != meshFaceI && layerFaceType[facei] == -1
                    && !boundaryFaces[facei] && cFaces.found(facei)
                )
                {
                    layerFaceType[facei] = 1;
                }
            }
        }

        forAll(f, fp)
        {
            label meshPointI = f[fp];
            const labelList& pFaces = mesh.pointFaces()[meshPointI];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if
                (
                    layerFaceType[facei] == -1
                    && !boundaryFaces[facei] && cFaces.found(facei)
                 )
                {
                    layerFaceType[facei] = 1;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        bPoints,
        orEqOp<bool>(),
        false
    );

    forAll(bPoints, meshPointI)
    {
        if (!bPoints[meshPointI])
        {
            continue;
        }

        const labelList& pCells = mesh.pointCells()[meshPointI];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];

            if (cellLayerNumber[celli] == 0)
            {
                cellLayerNumber[celli]++;
            }

            const labelList& pEdges = mesh.pointEdges()[meshPointI];

            labelHashSet cEdges(mesh.cellEdges()[celli]);
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if
                (
                    !boundaryEdges[edgei] && stackEdgeLayer[edgei] == 0
                    && cEdges.found(edgei)
                 )
                {
                    stackEdgeLayer[edgei]++;
                }
            }
            const labelList& pFaces = mesh.pointFaces()[meshPointI];
            labelHashSet cFaces(mesh.cells()[celli]);

            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if
                (
                    layerFaceType[facei] == -1
                    &&!boundaryFaces[facei] && cFaces.found(facei)
                 )
                {
                    layerFaceType[facei] = 1;
                }
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncEdgeList
        (
            mesh,
            stackEdgeLayer,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncFaceList
        (
            mesh,
            layerFaceType,
            maxEqOp<label>()
        );

        labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
        labelList neiCellLayerNumber(mesh.nFaces()-mesh.nInternalFaces());
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            neiLayerCells[facei-mesh.nInternalFaces()] =
                layerCells[owners[facei]];
            neiCellLayerNumber[facei-mesh.nInternalFaces()] =
                cellLayerNumber[owners[facei]];
        }
        syncTools::swapBoundaryFaceList(mesh, neiLayerCells);
        syncTools::swapBoundaryFaceList(mesh, neiCellLayerNumber);

        forAll(mesh.faces(), meshFaceI)
        {
            label own = owners[meshFaceI];

            label patchI = patches.whichPatch(meshFaceI);

            if (layerFaceType[meshFaceI] != -1)
            {
                continue;
            }

            if (patchI == -1 || patches[patchI].coupled())
            {
                label neiCLN = -1;
                label nLayerCell = -1;
                label nei = -1;

                if (patchI == -1)
                {
                    nei = mesh.faceNeighbour()[meshFaceI];
                    neiCLN = cellLayerNumber[nei];
                    nLayerCell = layerCells[nei];
                }
                else
                {
                    neiCLN =
                        neiCellLayerNumber[meshFaceI-mesh.nInternalFaces()];
                    nLayerCell =
                        neiLayerCells[meshFaceI-mesh.nInternalFaces()];
                }

                label celli = -1;
                label prevLayer = -1;

                if (layerCells[own] == nLayerCell)
                {
                    if (cellLayerNumber[own] > 0 && neiCLN == 0)
                    {
                        if (patchI == -1)
                        {
                            celli = nei;
                            prevLayer = cellLayerNumber[own];
                        }
                        else
                        {
                            nSet++;
                        }
                    }
                    else if
                    (
                        neiCLN > 0 && cellLayerNumber[own] == 0
                    )
                    {
                        celli = own;
                        prevLayer = neiCLN;
                    }
                }

                if (celli != -1)
                {
                    nSet++;
                    cellLayerNumber[celli] = prevLayer + 1;
                    layerFaceType[meshFaceI] = 2;

                    face f = mesh.faces()[meshFaceI];
                    const labelList& fEdges = mesh.faceEdges()[meshFaceI];
                    labelHashSet fEdgeSet(fEdges);

                    forAll(f, fp)
                    {
                        label meshPointI = f[fp];
                        const labelList& pEdges =
                            mesh.pointEdges()[meshPointI];

                        forAll(pEdges, pEI)
                        {
                            labelHashSet cEdges(mesh.cellEdges()[celli]);
                            forAll(pEdges, pEI)
                            {
                                label edgei = pEdges[pEI];
                                if
                                (
                                    stackEdgeLayer[edgei] == 0
                                    && cEdges.found(edgei)
                                    && !fEdgeSet.found(edgei)
                                )
                                {
                                    stackEdgeLayer[edgei]++;
                                }
                            }
                        }
                    }

                    labelHashSet cFaces(mesh.cells()[celli]);
                    forAll(fEdges, fEI)
                    {
                        label edgei = fEdges[fEI];

                        const labelList& eFaces = mesh.edgeFaces()[edgei];
                        forAll(eFaces, eFI)
                        {
                            label facei = eFaces[eFI];
                            if
                            (
                                facei != meshFaceI && layerFaceType[facei] == -1
                                && cFaces.found(facei)
                             )
                            {
                                layerFaceType[facei] = 1;
                            }
                        }
                    }
                    forAll(f, fp)
                    {
                        label meshPointI = f[fp];
                        const labelList& pFaces =
                            mesh.pointFaces()[meshPointI];
                        forAll(pFaces, pFI)
                        {
                            label facei = pFaces[pFI];
                            if
                            (
                                layerFaceType[facei] == -1
                                && cFaces.found(facei)
                            )
                            {
                                layerFaceType[facei] = 1;
                            }
                        }
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    //Mark outer stack
    forAll(mesh.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[facei];
            label nLayerCell = -1;

            if (patchI == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[facei-mesh.nInternalFaces()];
            }

            if
            (
                (layerCells[own] != nLayerCell)
                && layerFaceType[facei] == -1
            )
            {
                layerFaceType[facei] = 3;
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh,
        layerFaceType,
        maxEqOp<label>()
    );

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerPointType[meshPointI] = 1;
    }

    syncTools::syncPointList
    (
        mesh,
        layerPointType,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(mesh.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[facei];
            label nLayerCell = -1;

            if (patchI == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[facei-mesh.nInternalFaces()];
            }

            if (layerFaceType[facei] == 2)
            {
                const face f = mesh.faces()[facei];
                forAll(f,fp)
                {
                    if (layerPointType[f[fp]] == -1)
                    {
                        layerPointType[f[fp]] = 0;
                    }
                }
            }
            else if (layerFaceType[facei] == 3)
            {
                const face f = mesh.faces()[facei];
                if (layerCells[own] == -1 || nLayerCell == -1)
                {
                    forAll(f,fp)
                    {
                        layerPointType[f[fp]] = 2;
                    }
                }
            }
        }
    }

    forAll(mesh.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[facei];
            label nLayerCell = -1;
            if (patchI == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[facei-mesh.nInternalFaces()];
            }

            if (layerFaceType[facei] == 3)
            {
                const face f = mesh.faces()[facei];
                if (layerCells[own] != -1 && nLayerCell != -1)
                {
                    forAll(f,fp)
                    {
                        layerPointType[f[fp]] = 3;
                    }
                }
            }
        }
    }
    syncTools::syncPointList(mesh, layerPointType, maxEqOp<label>(), label(-1));

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            const labelList& cellPoints = mesh.cellPoints()[celli];
            forAll(cellPoints, cPtI)
            {
                label pointi = cellPoints[cPtI];

                if (layerPointType[pointi] == -1)
                {
                    layerPointType[cellPoints[cPtI]] = 3;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        layerPointType,
        maxEqOp<label>(),
        label(-1)
    );

    if (debug)
    {
        Time& runTime = const_cast<Time&>(mesh.time());
        runTime++;
        mesh.write();

        volScalarField layerCount
        (
            IOobject
            (
                "layerCount",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(layerCount, celli)
        {
            layerCount[celli] = cellLayerNumber[celli];
        }
        layerCount.write();
    }
}


void Foam::layerManipulate::writeLayerInfo()
{
    Info<<"Writing layer information "<<endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const labelIOList& layerCells = meshRefiner_.layerCells()();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    const labelList& meshPoints = pp.meshPoints();

    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    boolList boundaryFaces(mesh.nFaces(), false);
    boolList boundaryEdges(mesh.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    forAll(pp.edges(), edgei)
    {
        boundaryEdges[meshEdges[edgei]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh.nCells(), 0);
    labelList stackEdgeLayer(mesh.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh.nFaces(), -1);
    labelList layerPointType(mesh.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    boolList stackEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        if (stackEdgeLayer[edgei] != 0)
        {
            stackEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        pp,
        meshEdges,
        0.707,
        0.707
    );

    boolList excludedFaces(pp.size(), false);
    pointField pointNormals = eClass.calculatePointNormals
    (
        excludedFaces,
        0,
        true
    );

    labelList layerCount(mesh.nPoints(), -1);
    pointField pointOrigin(mesh.nPoints(), vector::zero);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = 0;
    }

    boolList edgeSet(mesh.nEdges(), false);
    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncEdgeList
        (
            mesh,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        forAll(mesh.edges(), edgei)
        {
            if (!isMasterEdge[edgei])
            {
                continue;
            }

            if (!boundaryEdges[edgei] && stackEdges[edgei] && !edgeSet[edgei])
            {
                edge e = mesh.edges()[edgei];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgei] = true;
                        nSet++;
                    }
                }

                if (!edgeSet[edgei] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgei] = true;
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    const labelList& owners = mesh.faceOwner();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList totalNumLayers(mesh.nPoints(), -1);
    forAll(mesh.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            label own = owners[facei];
            label layerCount = -1;

            if (layerCells[own] != -1)
            {
                layerCount = cellLayerNumber[own];
            }
            else
            {
                label patchI = patches.whichPatch(facei);
                if (patchI == -1)
                {
                    label nei = mesh.faceNeighbour()[facei];
                    if (layerCells[nei] != -1)
                    {
                        layerCount = cellLayerNumber[nei];
                    }
                }
            }

            if (layerCount != -1)
            {
                face f = mesh.faces()[facei];
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    totalNumLayers[pointi] = max
                    (
                        totalNumLayers[pointi],layerCount
                    );
                }
            }
        }
    }

    pointField outerPt(mesh.nPoints(), vector::zero);
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            outerPt[pointi] = mesh.points()[pointi];
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh,
            outerPt,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        forAll(mesh.edges(), edgei)
        {
            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    outerPt[e[1]] = outerPt[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    outerPt[e[0]] = outerPt[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    labelList patchNumLayers(pp.meshPoints().size(), 0);
    labelList patchLevel(pp.meshPoints().size(), 0);
    scalarField patchFCH(pp.meshPoints().size(), 0);
    scalarField patchLayerHeight(pp.meshPoints().size(), 0);

    forAll(pp.meshPoints(), pointi)
    {
        label meshPointI = pp.meshPoints()[pointi];
        patchNumLayers[pointi] = totalNumLayers[meshPointI];
        patchLevel[pointi] = pointLevel[meshPointI];

        const labelList& pEdges = mesh.pointEdges()[meshPointI];

        scalar firstCellHeight = GREAT;
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];

            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                vector eVec = mesh.edges()[edgei].vec(mesh.points());
                scalar projEdgeDist = mag(eVec & pointNormals[pointi]);
                firstCellHeight = min(firstCellHeight, projEdgeDist);
            }
        }
        patchFCH[pointi] = firstCellHeight;

        vector layerVec = outerPt[meshPointI] - mesh.points()[meshPointI];
        patchLayerHeight[pointi] = mag(layerVec & pointNormals[pointi]);
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        patchLayerHeight,
        minEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        patchFCH,
        minEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        patchLevel,
        maxEqOp<label>(),
        label(0)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        patchNumLayers,
        maxEqOp<label>(),
        label(0)          // null value
    );

    //Convert to face data for output
    labelList facePatchNumLayers(pp.size(), labelMax);
    labelList facePatchLevel(pp.size(), 0);
    scalarField facePatchFCH(pp.size(), 0);
    scalarField facePatchLayerHeight(pp.size(), 0);

    label nPatches = patches.size();
    scalarField avePatchNumLayer(nPatches, scalar(0));
    scalarField avePatchFCH(nPatches, scalar(0));
    scalarField avePatchThickness(nPatches, scalar(0));

    forAll(pp, facei)
    {
        label own = owners[pp.addressing()[facei]];
        if (layerCells[own] < 0)
        {
           continue;
        }

        const face& f = pp.localFaces()[facei];

        forAll(f, fp)
        {
            facePatchFCH[facei] += patchFCH[f[fp]];
            facePatchLayerHeight[facei] += patchLayerHeight[f[fp]];
            facePatchNumLayers[facei] = min
                (patchNumLayers[f[fp]],facePatchNumLayers[facei]);
            facePatchLevel[facei] = max
                (patchLevel[f[fp]],facePatchLevel[facei]);
        }
        facePatchFCH[facei] /= f.size();
        facePatchLayerHeight[facei] /= f.size();

        label patchi = patches.whichPatch(pp.addressing()[facei]);
        if (patchi != -1)
        {
            avePatchNumLayer[patchi] += facePatchNumLayers[facei];
            avePatchFCH[patchi] += facePatchFCH[facei];
            avePatchThickness[patchi] += facePatchLayerHeight[facei];
        }
    }

    boolList validPatches(patches.size(), false);
    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            label sz = patches[patchi].size();
            reduce(sz, sumOp<label>());
            if (sz > 0)
            {
                scalar avgNumLayers =
                    returnReduce(avePatchNumLayer[patchi], sumOp<scalar>())/sz;
                if (avgNumLayers > 0)
                {
                    validPatches[patchi] = true;
                }
            }
        }
    }

    // Find maximum length of a patch name, for a nicer output
    label maxPatchNameLen = 0;
    forAll(patches, patchi)
    {
        if (validPatches[patchi])
        {
            word patchName = patches[patchi].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }
    }

    Info<< nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
        << setw(0) << " faces    layers avg thickness[m]" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << " "
        << setw(0) << "                 near-wall overall" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
        << setw(0) << " -----    ------ --------- -------" << endl;

    forAll(patches, patchi)
    {
        if (validPatches[patchi])
        {
            label sz = patches[patchi].size();
            reduce(sz, sumOp<label>());

            if (sz > 0)
            {
                scalar avgThickness =
                    returnReduce(avePatchThickness[patchi], sumOp<scalar>())/sz;
                scalar avgNearWallThickness =
                    returnReduce(avePatchFCH[patchi], sumOp<scalar>())/sz;
                scalar avgNumLayers =
                    returnReduce(avePatchNumLayer[patchi], sumOp<scalar>())/sz;
                Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                    << patches[patchi].name() << setprecision(3)
                    << " " << setw(8) << sz
                    << " " << setw(6) << avgNumLayers
                    << " " << setw(8) << avgNearWallThickness
                    << "  " << setw(8) << avgThickness
                    << endl;
            }
        }
    }

    if (layerParams_.writeVTK() || layerParams_.writeExtraVTK())
    {
        simpleVTKWriter layerVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );
        layerVTK.addFaceData("numLayers", facePatchNumLayers);
        layerVTK.addFaceData("fch", facePatchFCH);
        layerVTK.addFaceData("level", facePatchLevel);
        layerVTK.addFaceData("layerHeight", facePatchLayerHeight);

        if (layerParams_.writeExtraVTK())
        {
            scalarField facePatchPID(pp.size(), -1);
            forAll(pp, i)
            {
                label patchi = patches.whichPatch(pp.addressing()[i]);
                facePatchPID[i] = patchi;
            }
            layerVTK.addFaceData("pid", facePatchPID);
        }

        fileName linfoPath("layerInfo"/mesh.name());
        if (Pstream::master())
        {
            if (!isDir(linfoPath))
            {
                mkDir(linfoPath);
            }
        }
        layerVTK.write(linfoPath/"layerInfo.vtk");
    }
}


void Foam::layerManipulate::setLayerMethod
(
    const indirectPrimitivePatch& pp,
    const labelList& meshPoints,
    const labelList& ppPtLevel,
    const labelList& ppPtLevelMin,
    List<Tuple2<word,scalar>>& layerMethod,
    scalarField& ppMaxLayerThickness
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const scalar edge0Length = meshRefiner_.meshCutter().level0EdgeLength();

    const List<wordList>& layerSpec = layerParams_.layerSpec();
    const scalarField& patchFCH = layerParams_.fch();
    const scalarField& patchExpansionRatio = layerParams_.expansionRatio();
    const scalarField& patchMaxLayerThickness =
        layerParams_.maxLayerThickness();
    const scalarField& patchFinalLayerThickness =
        layerParams_.finalLayerThickness();
    const labelList& patchNumLayers = layerParams_.numLayers();

    labelList syncType(meshPoints.size(), labelMax);
    scalarField syncFCH(meshPoints.size(), GREAT);
    scalarField syncExpansion(meshPoints.size(), GREAT);

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        label patchI = patches.whichPatch(meshFaceI);

        face f = pp.localFaces()[i];

        word method = layerSpec[patchI][1];
        if (method == "fch" || method == "rfch")
        {
            bool relativeFCH(method == "rfch" ? true : false);
            scalar relativeScaling = 1;
            forAll(f,fp)
            {
                if (relativeFCH)
                {
                    relativeScaling = edge0Length / (1<<ppPtLevel[f[fp]]);
                }
                label pointi = f[fp];
                syncType[pointi] = min(syncType[pointi],label(0));
                syncFCH[pointi] = min
                (
                    syncFCH[pointi],
                    relativeScaling*patchFCH[patchI]
                );
            }
        }
        else if (method == "expansionRatio")
        {
            forAll(f,fp)
            {
                label pointi = f[fp];
                syncType[pointi] = min(syncType[pointi],label(1));
                syncExpansion[pointi] = min
                (
                    syncExpansion[pointi],
                    patchExpansionRatio[patchI]
                );
            }
        }
        else
        {
            WarningInFunction
                << "Cannot find correct layer method for patch "
                << patches[patchI].name()
                << " requires fch, rfch or expansionRatio to be set"
                << endl;
        }

        word optionalConstraint = layerSpec[patchI][2];
        if
        (
            optionalConstraint == "maxLayerThicknessAbs"
            || optionalConstraint == "maxLayerThickness"
        )
        {
            bool relativeMLT = true;
            if (optionalConstraint == "maxLayerThicknessAbs")
            {
                relativeMLT = false;
            }
            scalar plt = patchMaxLayerThickness[patchI];
            if (relativeMLT)
            {
                plt *= edge0Length;
            }

            forAll(f,fp)
            {
                label pti = f[fp];
                scalar lLen = plt;
                if (relativeMLT)
                {
                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);

                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));
                    lLen *= levelWeight;
                }
                ppMaxLayerThickness[pti] = min
                (
                    ppMaxLayerThickness[pti],
                    lLen
                );
            }
        }
        else if
        (
            optionalConstraint == "finalLayerThicknessAbs"
            || optionalConstraint == "finalLayerThickness"
        )
        {
            bool relativeFLT = true;
            if (optionalConstraint == "finalLayerThicknessAbs")
            {
                relativeFLT = false;
            }
            if (method == "fch" || method == "rfch")
            {
                scalar pflt = patchFinalLayerThickness[patchI];
                scalar pfch = patchFCH[patchI];
                label pnl = patchNumLayers[patchI];

                bool relativeFCH(method == "rfch" ? true : false);
                forAll(f,fp)
                {
                    label pti = f[fp];
                    scalar ptFCH = pfch;
                    scalar ptFLT = pflt;

                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);
                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));

                    if (relativeFCH)
                    {
                        ptFCH *= edge0Length*levelWeight;
                    }
                    if (relativeFLT)
                    {
                        ptFLT *= edge0Length*levelWeight;
                    }

                    scalar str = 1.0;
                    if (pnl > 1)
                    {
                        str = pow((ptFLT/ptFCH),(1.0/(pnl-1)));
                    }
                    scalar lLen = pnl*ptFCH;
                    if (str < 1.0 - SMALL || str > 1.0 + SMALL)
                    {
                        lLen = ptFCH*(1.-pow(str,pnl))
                            / (1.-str);
                    }
                    ppMaxLayerThickness[pti] = min
                    (
                        ppMaxLayerThickness[pti],
                        lLen
                    );
                }
            }
            else
            {
                scalar pflt = patchFinalLayerThickness[patchI];
                label pnl = patchNumLayers[patchI];
                scalar invstr = 1.0/patchExpansionRatio[patchI];

                forAll(f,fp)
                {
                    label pti = f[fp];
                    scalar ptFLT = pflt;

                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);
                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));

                    if (relativeFLT)
                    {
                        ptFLT *= edge0Length*levelWeight;
                    }
                    scalar lLen = pnl*ptFLT;
                    if (invstr < 1.0 - SMALL || invstr > 1.0 + SMALL)
                    {
                        lLen = ptFLT*(1.-pow(invstr,pnl))
                            / (1.-invstr);
                    }
                    ppMaxLayerThickness[pti] = min
                    (
                        ppMaxLayerThickness[pti],
                        lLen
                    );
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        ppMaxLayerThickness,
        minEqOp<scalar>(),
        GREAT         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        syncType,
        minEqOp<label>(),
        label(-1)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        syncExpansion,
        minEqOp<scalar>(),
        GREAT          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        syncFCH,
        minEqOp<scalar>(),
        GREAT          // null value
    );

    forAll(meshPoints, pointi)
    {
       if (syncType[pointi] == 0)
       {
           layerMethod[pointi] = Tuple2<word,scalar>
               ("fch",syncFCH[pointi]);
       }
       else if (syncType[pointi] == 1)
       {
           layerMethod[pointi] = Tuple2<word,scalar>
               ("expansionRatio",syncExpansion[pointi]);
       }
    }

    return;
}


void Foam::layerManipulate::updateGeometry
(
    const pointField& newPoints,
    pointField& newFaceCentres,
    vectorField& newFaceAreas,
    pointField& newCellCentres,
    scalarField& newCellVolumes
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    newFaceCentres = vector::zero;
    newFaceAreas = vector::zero;
    newCellCentres = vector::zero;
    newCellVolumes = 0;

    if (layerParams_.fastGeomUpdate())
    {
        scalar a1 = (1.0/3.0);
        point avePt = Zero;
        label nextFp = 0;
        vector sumN = Zero;
        scalar sumA = 0.0;
        vector sumAc = Zero;

        vector n;
        vector c;
        scalar a;

        const faceList& fs = mesh.faces();

        forAll(fs, facei)
        {
            const face& f = fs[facei];
            label nPoints = f.size();
            point pt0 = newPoints[f[0]];
            pointField facePts(nPoints,vector::zero);
            forAll(f,fp)
            {
               facePts[fp] = newPoints[f[fp]]-pt0;
            }
            if (nPoints == 3)
            {
                newFaceCentres[facei] = (1.0/3.0)*
                (
                    facePts[0] + facePts[1] + facePts[2]
                );
                newFaceAreas[facei] = 0.5*
                (
                    (facePts[1] - facePts[0]) ^(facePts[2] - facePts[0])
                );
            }
            else
            {
                sumN = Zero;
                sumA = 0.0;
                sumAc = Zero;
                const point& pt0 = facePts[0];
                nextFp = 1;
                for (label fp = 0; fp < nPoints-2; fp++)
                {
                    const point& pt1 = facePts[nextFp];
                    nextFp = f.fcIndex(nextFp);
                    const point& pt2 = facePts[nextFp];
                    n = (pt1 - pt0)^(pt2 -pt0);
                    c = pt0 + pt1 +pt2;
                    a = mag(n);
                    sumA += a;
                    sumAc += a*c;
                    sumN += n;
                }
                if (sumA > VSMALL)
                {
                    newFaceCentres[facei] = a1*sumAc/sumA;
                }
                else
                {
                    avePt = Zero;
                    forAll(f,fp)
                    {
                        avePt += facePts[fp];
                    }
                    newFaceCentres[facei] = (avePt/nPoints);
                }
                newFaceAreas[facei] = 0.5*sumN;
            }
            newFaceCentres[facei] += pt0;
        }

        const labelList& owners = mesh.faceOwner();
        const labelList& neighbours = mesh.faceNeighbour();
        pointField cEst(mesh.nCells(), Zero);
        forAll(owners, facei)
        {
            label own = owners[facei];
            cEst[own] += newFaceCentres[facei];
        }
        forAll(neighbours, facei)
        {
            label nei = neighbours[facei];
            cEst[nei] += newFaceCentres[facei];
        }
        forAll(cEst, celli)
        {
            cEst[celli] /= mesh.cells()[celli].size();
        }

        scalarField cellCentreMag(mesh.nCells(), 0);
        forAll(owners, facei)
        {
            label own = owners[facei];
            const point& cc = cEst[own];
            const point& fC = newFaceCentres[facei];
            const point& fA = newFaceAreas[facei];
            scalar  pyr3Vol = fA & (fC - cc);
            vector pc = 0.75*fC + 0.25*cc;
            scalar pyr3VolMag = mag(pyr3Vol);
            newCellCentres[own] += pyr3VolMag*pc;
            // Accumulate face-pyramid volume
            newCellVolumes[own] += pyr3Vol;
            cellCentreMag[own] += pyr3VolMag;
         }

        forAll(neighbours, facei)
        {
            label nei = neighbours[facei];
            const point& cc = cEst[nei];
            const point& fC = newFaceCentres[facei];
            const point& fA = newFaceAreas[facei];
            scalar  pyr3Vol = fA & (cc - fC);
            vector pc = 0.75*fC + 0.25*cc;
            scalar pyr3VolMag = mag(pyr3Vol);
            newCellCentres[nei] += pyr3VolMag*pc;
            // Accumulate face-pyramid volume
            newCellVolumes[nei] += pyr3Vol;
            cellCentreMag[nei] += pyr3VolMag;
        }

        forAll(newCellCentres, celli)
        {
            scalar cellVol = cellCentreMag[celli];
            if (cellVol > VSMALL)
            {
                newCellCentres[celli] /= cellVol;
            }
            else
            {
                newCellCentres[celli] = cEst[celli];
            }
        }

        newCellVolumes  *= (1.0/3.0);
    }
    else
    {
        scalarField newMagFaceAreas(mesh.nFaces(), Zero);

        // Geometrical calculations
        mesh.makeFaceCentresAndAreas
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newMagFaceAreas
        );

        mesh.makeCellCentresAndVols
        (
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );
    }

    return;
}


void Foam::layerManipulate::deactivateSquishCells
(
    const labelList& layerPointType,
    const labelList& boundaryPoints,
    const pointField& newFaceCentres,
    const vectorField& newFaceAreas,
    const pointField& newCellCentres,
    labelList& stationaryPts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelIOList& layerCells = meshRefiner_.layerCells()();
    const scalar& squishTol = layerParams_.squishTol();

    forAll(stationaryPts, pointi)
    {
        if (stationaryPts[pointi] == 5)
        {
            stationaryPts[pointi] = -2;
        }
    }

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            continue;
        }

        label nLayerPts = 0;
        label nZonePts = 0;
        label nBoundaryPts = 0;
        const labelList& cPts = mesh.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];
            if (layerPointType[pointi] != -1)
            {
                nLayerPts++;
            }
            //zone pts
            if
            (
                stationaryPts[pointi] == -1
                || stationaryPts[pointi] == 1
            )
            {
                nZonePts++;
            }
            //boundary pts
            if (boundaryPoints[pointi] > -1)
            {
                nBoundaryPts++;
            }
        }
        if (nLayerPts > 0 || nZonePts > 0 || nBoundaryPts > 0)
        {
            continue;
        }

        bool squished = false;
        const cell& c = mesh.cells()[celli];
        scalar squishVal = scalar(0);
        scalar totFaceArea = scalar(0);
        forAll(c, cFI)
        {
            label facei = c[cFI];
            vector fcToCC = newFaceCentres[facei] - newCellCentres[celli];
            scalar magfcToCC = mag(fcToCC);
            vector faceArea = newFaceAreas[facei];
            scalar magFaceArea = mag(faceArea);
            totFaceArea += magFaceArea;
            if (magFaceArea < VSMALL || magfcToCC < VSMALL)
            {
                squished = true;
            }
            else
            {
                fcToCC /= magfcToCC;
                squishVal += mag(faceArea & fcToCC);
            }
        }
        if (totFaceArea > VSMALL)
        {
            squishVal /= totFaceArea;
            if (squished || squishVal < squishTol)
            {
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    if (stationaryPts[pointi] == -2)
                    {
                        stationaryPts[pointi] = 5;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    DynamicList<label> nbrCells(mesh.nCells()/10);
    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            continue;
        }

        label nLayerPts = 0;
        label nStaticPts = 0;
        const labelList& cPts = mesh.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];

            if (stationaryPts[pointi] == 5)
            {
                nStaticPts++;
            }
            if (layerPointType[pointi] != -1)
            {
                nLayerPts++;
            }
        }
        if (nLayerPts > 0 && nStaticPts > 0)
        {
            nbrCells.append(celli);
        }
    }

    forAll(nbrCells, nbri)
    {
        label celli = nbrCells[nbri];
        const labelList& cPts = mesh.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];
            if (stationaryPts[pointi] == -2)
            {
                stationaryPts[pointi] = 5;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );
}

void Foam::layerManipulate::smoothLayerStack
(
    const meshControl& controller,
    const label nSmoothIter,
    const dictionary& grownUpGeometryDict,
    const dictionary& grownUpZoneGeometryDict,
    const labelList& layerOffset,
    labelList& adjustedNLayers,
    const bool adjustFinalNLayers
)
{
    Info<<"Smoothing layer cells "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const scalar edge0Length = meshRefiner_.meshCutter().level0EdgeLength();

    const labelIOList& layerCells = meshRefiner_.layerCells()();

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    boolList layerEdges(mesh.nEdges(), false);
    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            const labelList& cellEdges = mesh.cellEdges()[celli];
            forAll(cellEdges, cEI)
            {
                layerEdges[cellEdges[cEI]] = true;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    //Set layer growth method
    labelList ppPtLevel(meshPoints.size(), -1);
    labelList ppPtLevelMin(meshPoints.size(), labelMax);
    forAll(meshPoints, pti)
    {
        const labelList& pFaces = pp.pointFaces()[pti];
        forAll(pFaces, pfi)
        {
            label meshFaceI = pp.addressing()[pFaces[pfi]];
            label own = mesh.faceOwner()[meshFaceI];
            ppPtLevel[pti] = max(cellLevel[own],ppPtLevel[pti]);
            ppPtLevelMin[pti] = min(cellLevel[own],ppPtLevelMin[pti]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        ppPtLevel,
        maxEqOp<label>(),
        label(-1)         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        ppPtLevelMin,
        minEqOp<label>(),
        labelMax         // null value
    );

    List<Tuple2<word,scalar>> layerMethod(meshPoints.size());
    scalarField ppMaxLayerThickness(meshPoints.size(), GREAT);
    setLayerMethod
    (
        pp,
        meshPoints,
        ppPtLevel,
        ppPtLevelMin,
        layerMethod,
        ppMaxLayerThickness
    );

    boolList boundaryFaces(mesh.nFaces(), false);
    boolList boundaryEdges(mesh.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    autoPtr<indirectPrimitivePatch> grownUpSnapPPPtr = makeGrownUpPatch(true);
    indirectPrimitivePatch& grownUpSnapPP = grownUpSnapPPPtr();

    autoPtr<searchableSurfaces> grownUpGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                mesh.time().constant(),     // directory
                "triSurface",               // instance
                mesh.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
             ),
             grownUpGeometryDict
         )
    );

    autoPtr<searchableSurfaces> grownUpZoneGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                mesh.time().constant(),     // directory
                "triSurface",               // instance
                mesh.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
             ),
             grownUpZoneGeometryDict
         )
    );

    autoPtr<indirectPrimitivePatch> grownUpPPPtr = makeGrownUpPatch();
    indirectPrimitivePatch& grownUpPP = grownUpPPPtr();

    autoPtr<indirectPrimitivePatch> fzonePPPtr = makeFaceZonePatch();
    indirectPrimitivePatch& fzonePP = fzonePPPtr();

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
    boolList featurePts(mesh.nPoints(), false);
    boolList featureEdges(mesh.nEdges(), false);

    //Mark zone feature edges
    if
    (
        layerParams_.dualReSnapZones()
        && returnReduce(fzonePP.size(), sumOp<label>()) != 0
    )
    {
        const labelList zoneEdges
        (
            fzonePP.meshEdges
            (
                mesh.edges(),
                mesh.pointEdges()
            )
        );

        autoPtr<boolList> flipMap
        (
            new boolList(fzonePP.size(), false)
        );
        labelList neiCellZone;
        labelList cellToZone(mesh.nCells(), -1);
        forAll(mesh.cells(), celli)
        {
            label zonei = mesh.cellZones().whichZone(celli);
            cellToZone[celli] = zonei;
        }
        syncTools::swapBoundaryCellList(mesh, cellToZone, neiCellZone);

        forAll(fzonePP, i)
        {
            label facei = fzonePP.addressing()[i];
            label ownZoneI = mesh.cellZones().whichZone
            (
                mesh.faceOwner()[facei]
            );

            label  neiZone = -1;
            if (mesh.isInternalFace(facei))
            {
                neiZone = mesh.cellZones().whichZone
                (
                    mesh.faceNeighbour()[facei]
                );
            }
            else
            {
                neiZone = neiCellZone[facei-mesh.nInternalFaces()];
            }
            if (ownZoneI < neiZone)
            {
                flipMap()[i] = true;
            }
        }

        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            fzonePP,
            zoneEdges,
            0.8191,
            0.8191,
            -GREAT,
            flipMap // optional to include a flip map
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();
        forAll(zoneEdges, edgei)
        {
            if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
                || eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = zoneEdges[edgei];
                edge e = mesh.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }
    }

    const labelList meshEdgesGrownUp
    (
        grownUpPP.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    edgeClassification eClassGrownUp
    (
        mesh,
        mesh.points(),
        grownUpPP,
        meshEdgesGrownUp,
        0.8191,
        0.8191
    );

    autoPtr<indexedOctree<treeDataEdge>> grownUpFeatureMeshPtr;
    edgeMesh eMeshGrownUp =
       meshRefiner_.calcFeatureEdgeMesh(mesh,grownUpPP,eClassGrownUp);

    // Calculate grown up patch mesh fetaure edges
    if (eMeshGrownUp.points().size() > 0)
    {
        Info<<"Constructing grown up patch feature edges " << nl <<endl;
        // Calculate bb of all points
        treeBoundBox bb(eMeshGrownUp.points());
        bb = bb.extend(1e-4);

        grownUpFeatureMeshPtr.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,
                    eMeshGrownUp.edges(),
                    eMeshGrownUp.points(),
                    identity(eMeshGrownUp.edges().size())
                ),
                bb,     // overall search domain
                8,      // maxLevel
                10,     // leafsize
                3.0     // du
            )
        );
    }

    boolList flushGrownUpPatch(grownUpPP.meshPoints().size(), false);
    //Calculate feature edges grown up patch
    {
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClassGrownUp.edgeTypes();

        forAll(meshEdgesGrownUp, edgei)
        {
            if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
                || eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = meshEdgesGrownUp[edgei];
                edge e = mesh.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            featureEdges,
            orEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh,
            featurePts,
            orEqOp<bool>(),
            false
        );

        labelList nFeatureEdges(mesh.nPoints(), 0);
        forAll(mesh.points(), pointi)
        {
            labelList pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                if (isMasterEdge[pEdges[pEI]] && featureEdges[pEdges[pEI]])
                {
                    nFeatureEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nFeatureEdges,
            plusEqOp<label>(),
            label(0)
        );

        forAll(mesh.points(), pointi)
        {
            if (boundaryPoints[pointi] == -1)
            {
                if (featurePts[pointi])
                {
                    if (nFeatureEdges[pointi] > 2)
                    {
                        boundaryPoints[pointi] = 3;
                    }
                    else
                    {
                        boundaryPoints[pointi] = 2;
                    }
                }
            }
        }

        //Check for flush patches to constrain
        if (layerParams_.flushGrownUpIDs().size() > 0)
        {
            const List<labelList> gUpEdgePatches =
                eClassGrownUp.uniqueEdgePatches();
            labelHashSet flushPatchIDs(layerParams_.flushGrownUpIDs());
            forAll(gUpEdgePatches, edgei)
            {
                if
                (
                    gUpEdgePatches[edgei].size() > 1
                    && eType[edgei].first() != edgeClassification::CONCAVE
                    && eType[edgei].first() != edgeClassification::CONVEX
                )
                {
                    forAll(gUpEdgePatches[edgei], gupid)
                    {
                        if (flushPatchIDs.found(gUpEdgePatches[edgei][gupid]))
                        {
                            const edge e = grownUpPP.edges()[edgei];
                            flushGrownUpPatch[e[0]] = true;
                            flushGrownUpPatch[e[1]] = true;
                            break;
                        }
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                grownUpPP.meshPoints(),
                flushGrownUpPatch,
                orEqOp<bool>(),
                false
            );
        }

        forAll(grownUpPP.meshPoints(), ptI)
        {
            label pointi = grownUpPP.meshPoints()[ptI];

            if (boundaryPoints[pointi] == -1)
            {
                boundaryPoints[pointi] = 1;
            }
        }
    }

    forAll(pp.edges(), edgei)
    {
        boundaryEdges[meshEdges[edgei]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh.nCells(), 0);
    labelList stackEdgeLayer(mesh.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh.nFaces(), -1);
    labelList layerPointType(mesh.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    autoPtr<indirectPrimitivePatch> outerShellPtr;
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        DynamicList<label> extFaces(mesh.nFaces()/10);
        forAll(mesh.faces(), facei)
        {
            if (layerFaceType[facei] == 3 && isMasterFace[facei])
            {
                extFaces.append(facei);
            }
        }
        outerShellPtr.reset
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), extFaces),
                mesh.points()
             )
        );
    }

    //Calculate layer interface weights based on number of layer
    // and non-layers cells in order to reduce jagged layer interface
    scalarField nNonLayerCells(mesh.nPoints(),0);
    scalarField nLayerCells(mesh.nPoints(),0);
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            if (boundaryPoints[pointi] == 1)
            {
                const labelList& pFaces = mesh.pointFaces()[pointi];
                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];
                    label patchI = patches.whichPatch(facei);

                    if (patchI != -1 && !patches[patchI].coupled())
                    {
                        label celli = mesh.faceOwner()[facei];
                        if (layerCells[celli] > -1)
                        {
                            nLayerCells[pointi] += scalar(1);
                        }
                        else
                        {
                            nNonLayerCells[pointi] += scalar(1);
                        }
                    }
                }
            }
            else
            {
                const labelList& pCells = mesh.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    if (layerCells[celli] > -1)
                    {
                        nLayerCells[pointi] += scalar(1);
                    }
                    else
                    {
                        nNonLayerCells[pointi] += scalar(1);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nLayerCells,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nNonLayerCells,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    boolList stackEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        if (stackEdgeLayer[edgei] != 0)
        {
            stackEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    pointField pointNormals,unsmoothedPointNormals;
    labelList ftrPointOrigin(mesh.nPoints(), -1);
    //Mark feature points
    {
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            0.707,//convex
            0.93969,//concave
            -0.707//-0.9848 //baffle
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(meshEdges, edgei)
        {
            edge e = mesh.edges()[meshEdges[edgei]];
            if (eType[edgei].first() == edgeClassification::BAFFLE)
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(2));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(2));
            }
            else if
            (
                eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(1));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(1));
            }
            else if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
            )
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(0));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(0));
            }
        }
    }

    //Calculate area weighted point normals (feature points fixed)
    {
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            0.707,//-0.5,
            0.707//-0.5
        );

        boolList excludedFaces(pp.size(), false);
        pointNormals = eClass.calculatePointNormals
        (
            excludedFaces,
            layerParams_.nSmoothNormalsExtrude(),
            true
        );
        unsmoothedPointNormals =
            eClass.calculatePointNormals(excludedFaces, 0, true);

    }

    syncTools::syncPointList
    (
        mesh,
        ftrPointOrigin,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(ftrPointOrigin, pointi)
    {
        if (ftrPointOrigin[pointi] == 2)
        {
            const labelList& pFaces = mesh.pointFaces()[pointi];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if (!boundaryFaces[facei] && mesh.faces()[facei].size() == 3)
                {
                    ftrPointOrigin[pointi] = -1;
                    break;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        ftrPointOrigin,
        minEqOp<label>(),
        label(-1)
    );

    mesh.clearOut();
    label nMedialAxisIter = mesh.globalData().nTotalPoints();

    // Distance to wall
    List<pointData> pointWallDist(mesh.nPoints());
    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    // 1. Calculate distance to points where displacement is specified.
    {
        labelList nSurfLayer(mesh.nPoints(), -1);
        forAll(mesh.points(), pointi)
        {
            if (layerPointType[pointi] >= 2)
            {
                nSurfLayer[pointi] = 0;
            }
        }

        while (true)
        {
            label nSet = 0;

            syncTools::syncPointList
            (
                mesh,
                nSurfLayer,
                maxEqOp<label>(),
                label(-1)              // null value
            );

            forAll(mesh.edges(), edgei)
            {
                if (!boundaryEdges[edgei] && stackEdges[edgei])
                {
                    edge e = mesh.edges()[edgei];

                    if
                    (
                        nSurfLayer[e[0]] != -1 && nSurfLayer[e[1]] == -1
                        && boundaryPoints[e[0]] != 0
                    )
                    {
                        nSurfLayer[e[1]] = nSurfLayer[e[0]]+1;
                        nSet++;
                    }
                    else if
                    (
                        nSurfLayer[e[1]] != -1 && nSurfLayer[e[0]] == -1
                        && boundaryPoints[e[1]] != 0
                    )
                    {
                        nSurfLayer[e[0]] = nSurfLayer[e[1]]+1;
                        nSet++;
                    }
                }
            }

            if (returnReduce(nSet, sumOp<label>()) == 0)
            {
                break;
            }
        }

        scalarField estimatedLayerHeight(meshPoints.size(), 0);

        scalarField len(meshPoints.size(), 0);
        labelList nPointFaces(meshPoints.size(), 0);

        forAll(meshPoints, ptI)
        {
            const labelList& pFaces = pp.pointFaces()[ptI];

            forAll(pFaces, pFI)
            {
                label meshFaceI = pp.addressing()[pFaces[pFI]];
                label own = mesh.faceOwner()[meshFaceI];
                len[ptI] += edge0Length / (1<<cellLevel[own]);
            }
            nPointFaces[ptI] += pFaces.size();
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            len,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nPointFaces,
            plusEqOp<label>(),
            label(0)          // null value
        );

        //assume stretching for total layer thickness estimate
        // needs to be made exact
        scalar str = 0.769;//0.869;
        forAll(meshPoints, ptI)
        {
            scalar ratio = (1-pow(str, nSurfLayer[meshPoints[ptI]]))
                / (scalar(1.)-str);
            estimatedLayerHeight[ptI] = ratio*(len[ptI]/nPointFaces[ptI]);
        }

        // Seed data.
        DynamicList<pointData> wallInfo(meshPoints.size());
        DynamicList<label> wallPoints(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];
            wallPoints.append(meshPointI);
            wallInfo.append
            (
                pointData
                (
                    mesh.points()[meshPointI],
                    0.0,
                    estimatedLayerHeight[patchPointI],
                    unsmoothedPointNormals[patchPointI]
                )
            );
        }
        wallInfo.shrink();
        wallPoints.shrink();

        List<pointData> edgeWallDist(mesh.nEdges());
        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh,
            wallPoints,
            wallInfo,

            pointWallDist,
            edgeWallDist,
            0,
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);
    }

    scalarField medialDist(mesh.nPoints(), scalar(0));
    scalarField medialLayerHeight(mesh.nPoints(), scalar(0));

    // 1. Medial axis points
    {
        const edgeList& edges = mesh.edges();
        List<pointData> pointMedialDist(mesh.nPoints());
        List<pointData> edgeMedialDist(mesh.nEdges());

        // Seed point data.
        DynamicList<pointData> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        forAll(edges, edgei)
        {
            if (isMaxEdge(mesh, pointWallDist, edgei, -0.5))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                //const point medialAxisPt = e.centre(points);
                const edge& e = edges[edgei];
                vector eVec = e.vec(mesh.points());
                scalar eMag = mag(eVec);
                if (eMag > VSMALL)
                {
                    eVec /= eMag;

                    // Calculate distance along edge
                    const point& p0 = mesh.points()[e[0]];
                    const point& p1 = mesh.points()[e[1]];
                    scalar dist0 = (p0-pointWallDist[e[0]].origin()) & eVec;
                    scalar dist1 = (pointWallDist[e[1]].origin()-p1) & eVec;
                    scalar s = 0.5*(dist1+eMag+dist0);

                    point medialAxisPt;
                    if (s <= dist0)
                    {
                        medialAxisPt = p0;
                    }
                    else if (s >= dist0+eMag)
                    {
                        medialAxisPt = p1;
                    }
                    else
                    {
                        medialAxisPt = p0+(s-dist0)*eVec;
                    }

                    scalar totalLayerHeight = pointWallDist[e[0]].s() +
                        pointWallDist[e[1]].s();
                    vector origToOrig = pointWallDist[e[1]].origin()
                        - pointWallDist[e[0]].origin();

                    forAll(e, ep)
                    {
                        label pointi = e[ep];
                        maxPoints.append(pointi);
                        maxInfo.append
                        (
                            pointData
                            (
                                medialAxisPt,   //points[pointi],
                                magSqr(mesh.points()[pointi]-medialAxisPt),
                                totalLayerHeight,
                                origToOrig
                            )
                        );
                    }
                }
            }
        }
        maxInfo.shrink();
        maxPoints.shrink();

        // Do all calculations
        PointEdgeWave<pointData> medialDistCalc
        (
            mesh,
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            0,
            dummyTrackData
        );
        medialDistCalc.iterate(2*nMedialAxisIter);

        // Extract medial axis distance as pointScalarField
        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            if (pointMedialDist[meshPointI].valid(dummyTrackData))
            {
                scalar pLen = 2*edge0Length / (1<<ppPtLevel[ptI]);
                scalar maxDist = max(mag(pointMedialDist[meshPointI].v()),pLen);
                scalar mDist =
                    Foam::sqrt(pointMedialDist[meshPointI].distSqr());
                if (mDist < maxDist)
                {
                    medialDist[meshPointI] = mDist;
                }
                else
                {
                    medialDist[meshPointI] = GREAT;
                }
            }
            else
            {
                medialDist[meshPointI] = GREAT;
            }
        }

        forAll(mesh.points(), pointi)
        {
            if (pointMedialDist[pointi].valid(dummyTrackData))
            {
                medialLayerHeight[pointi] = pointMedialDist[pointi].s();
            }
            else
            {
                medialLayerHeight[pointi] = GREAT;
            }
        }
    }

    if (debug)
    {
        simpleVTKWriter medialVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );
        scalarField medialHeight(meshPoints.size(), 0);

        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            medialHeight[ptI] = medialDist[meshPointI];
        }
        medialVTK.addPointData("medialHeight", medialHeight);

        medialVTK.write("medialVTK.vtk");
    }

    labelList layerCount(mesh.nPoints(), -1);

    pointField pointOrigin(mesh.nPoints(), vector::zero);
    pointField pointSurfNormal(mesh.nPoints(), vector::zero);

    scalarField maxLayerThickness(mesh.nPoints(), GREAT);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = layerOffset[ptI];
        pointOrigin[meshPointI] = pp.localPoints()[ptI];
        pointSurfNormal[meshPointI] = -pointNormals[ptI];
        maxLayerThickness[meshPointI] = ppMaxLayerThickness[ptI];
    }

    labelList nVisited(mesh.nPoints(), 0);

    boolList edgeSet(mesh.nEdges(), false);

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            ftrPointOrigin,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh,
            pointOrigin,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointSurfNormal,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh,
            medialDist,
            maxEqOp<scalar>(),
            scalar(0)
        );

        syncTools::syncEdgeList
        (
            mesh,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        forAll(mesh.edges(), edgei)
        {
            if (!isMasterEdge[edgei])
            {
                continue;
            }

            if (!boundaryEdges[edgei] && stackEdges[edgei] && !edgeSet[edgei])
            {
                edge e = mesh.edges()[edgei];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgei] = true;

                        nVisited[e[1]]++;

                        maxLayerThickness[e[1]] = maxLayerThickness[e[0]];
                        ftrPointOrigin[e[1]] = ftrPointOrigin[e[0]];
                        pointOrigin[e[1]] += pointOrigin[e[0]];
                        medialDist[e[1]] += medialDist[e[0]];
                        pointSurfNormal[e[1]] += pointSurfNormal[e[0]];
                        nSet++;
                    }
                }

                if (!edgeSet[edgei] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgei] = true;

                        nVisited[e[0]]++;

                        maxLayerThickness[e[0]] = maxLayerThickness[e[1]];
                        ftrPointOrigin[e[0]] = ftrPointOrigin[e[1]];
                        pointOrigin[e[0]] += pointOrigin[e[1]];
                        medialDist[e[0]] += medialDist[e[1]];
                        pointSurfNormal[e[0]] += pointSurfNormal[e[1]];
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nVisited,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.points(), pointi)
    {
        if (nVisited[pointi] > 1)
        {
            pointOrigin[pointi] /= nVisited[pointi];
            pointSurfNormal[pointi] /= nVisited[pointi];
            medialDist[pointi] /= nVisited[pointi];
        }
    }

    //Check for unitinitilised hanging points
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh.points(), pointi)
        {
            if (layerPointType[pointi] == 2 && nVisited[pointi] == 0)
            {
                const labelList& pEdges = mesh.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    edge e = mesh.edges()[edgei];
                    label otherPt(e[0] ==  pointi ? e[1] : e[0]);
                    if (layerEdges[edgei] && layerPointType[otherPt] == 2)
                    {
                        ftrPointOrigin[pointi] =
                            max(ftrPointOrigin[otherPt],ftrPointOrigin[pointi]);
                        maxLayerThickness[pointi] = min
                        (
                            maxLayerThickness[otherPt],
                            maxLayerThickness[pointi]
                        );
                        pointOrigin[pointi] += pointOrigin[otherPt];
                        pointSurfNormal[pointi] += pointSurfNormal[otherPt];
                        medialDist[pointi] += medialDist[otherPt];
                        nVisited[pointi]++;
                    }
                }
            }
            else
            {
                nVisited[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nVisited,
            plusEqOp<label>(),
            label(0)
        );
        syncTools::syncPointList
        (
            mesh,
            ftrPointOrigin,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        syncTools::syncPointList
        (
            mesh,
            pointOrigin,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointSurfNormal,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            medialDist,
            plusEqOp<scalar>(),
            scalar(0)
        );

        forAll(mesh.points(), pointi)
        {
            if (nVisited[pointi] > 1)
            {
                pointOrigin[pointi] /= nVisited[pointi];
                pointSurfNormal[pointi] /= nVisited[pointi];
                medialDist[pointi] /= nVisited[pointi];
            }
        }
    }

    labelList totalNumLayers(mesh.nPoints(), -1);
    forAll(mesh.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            face f = mesh.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                totalNumLayers[pointi] = max
                (
                    totalNumLayers[pointi],layerCount[pointi]
                );
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(mesh.edges(), edgei)
        {
            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    pointField newPoints = mesh.points();
    pointField tPts(mesh.points());

    pointField newCellCentres(mesh.nCells(), vector::zero);
    scalarField newCellVolumes(mesh.nCells(), 0);
    pointField newFaceCentres(mesh.nFaces(), vector::zero);
    vectorField newFaceAreas(mesh.nFaces(), vector::zero);

    labelList resetPts(mesh.nPoints(), -1);

    //Prevent faceZones from moving if further than estimated
    // layer depth from surface
    // startionaryPts: -2 (moving) -1 (moving zone),
    // 1 (stationary), 2 (gap), 3 (non layer medial pts)
    labelList stationaryPts(mesh.nPoints(), -2);
    scalar zoneLayersScaling
    (
        layerParams_.dualReSnapZones()
        ? GREAT : layerParams_.zoneLayersScaling()
    );
    forAll(fzonePP.meshPoints(), pointi)
    {
        label meshPointI = fzonePP.meshPoints()[pointi];
        if (!pointWallDist[meshPointI].valid(dummyTrackData))
        {
            stationaryPts[meshPointI] = -1;
            continue;
        }

        scalar maxLayerHeightSqr = zoneLayersScaling
            *sqr(pointWallDist[meshPointI].s());
        if
        (
            pointWallDist[meshPointI].distSqr() == GREAT
            || pointWallDist[meshPointI].distSqr() >= maxLayerHeightSqr
        )
        {
            stationaryPts[meshPointI] = 1;
        }
        else
        {
            stationaryPts[meshPointI] = -1;
        }
    }

    forAll(grownUpPP.meshPoints(), pointi)
    {
        label meshPointI = grownUpPP.meshPoints()[pointi];
        if
        (
           !pointWallDist[meshPointI].valid(dummyTrackData)
           || !flushGrownUpPatch[pointi]
        )
        {
            continue;
        }

        scalar maxLayerHeightSqr = sqr(pointWallDist[meshPointI].s());
        if
        (
            pointWallDist[meshPointI].distSqr() == GREAT
            || pointWallDist[meshPointI].distSqr() >= maxLayerHeightSqr
        )
        {
            stationaryPts[meshPointI] = 1;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    //Check for cells with all points zoned and set stationary
    forAll(mesh.cells(), celli)
    {
        const labelList& cPts = mesh.cellPoints()[celli];
        label nZonePts = 0;
        forAll(cPts, i)
        {
            label meshPointI = cPts[i];
            if
            (
                stationaryPts[meshPointI] == -1
                || stationaryPts[meshPointI] == 1
            )
            {
                nZonePts++;
            }
        }
        if (nZonePts == cPts.size())
        {
            forAll(cPts, i)
            {
                label meshPointI = cPts[i];
                stationaryPts[meshPointI] = 1;
            }
        }
    }

    //Block based on medial axis edge cells
    {
        forAll(mesh.edges(), edgei)
        {
            if (isMaxEdge(mesh, pointWallDist, edgei, -0.5))
            {
                const labelList& eCells = mesh.edgeCells()[edgei];
                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    if (layerCells[celli] > -1)
                    {
                        continue;
                    }

                    const labelList& cPts = mesh.cellPoints()[celli];
                    label nLayerPts = 0;

                    forAll(cPts, i)
                    {
                        label pointi = cPts[i];
                        if (layerPointType[pointi] == -1)
                        {
                            nLayerPts++;
                        }
                    }

                    if (nLayerPts == cPts.size())
                    {
                        forAll(cPts, i)
                        {
                            label pointi = cPts[i];
                            if (stationaryPts[pointi] != 1)
                            {
                                stationaryPts[pointi] = 3;
                            }
                        }
                    }
                }
            }
        }
    }

    //Prevent motion of gap cell points
    if (mesh.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh.lookupObject<volScalarField>("gapCells");

        forAll(mesh.cells(), celli)
        {
            if (gapCells[celli] > -1)
            {
                const labelList& cPts = mesh.cellPoints()[celli];
                forAll(cPts, i)
                {
                    stationaryPts[cPts[i]] = 2;
                }
            }
        }
    }

    //Prevent motion at points further than the maximum projection distance
    const scalar maxProjDist = layerParams_.maxProjectionDist();
    if (maxProjDist < GREAT)
    {
        forAll(stationaryPts,pointi)
        {
            if
            (
                pointWallDist[pointi].valid(dummyTrackData)
                && Foam::sqrt(pointWallDist[pointi].distSqr()) > maxProjDist
            )
            {
                stationaryPts[pointi] = 4;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    boolList hangingNodes(mesh.nPoints(), false);
    List<labelList> hangingNbrs(mesh.nPoints(), labelList());
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh.cells(), celli)
        {
            const labelList& cEdges = mesh.cellEdges()[celli];
            const labelList& cPts = mesh.cellPoints()[celli];
            labelHashSet cEdgesSet(cEdges);
            label cLevel = cellLevel[celli];

            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                label pLevel = pointLevel[pointi];
                if (pLevel > cLevel)
                {
                    const labelList& pEdges = mesh.pointEdges()[pointi];
                    DynamicList<label> oPts(pEdges.size());
                    DynamicList<label> hPts(pEdges.size());

                    if (boundaryPoints[pointi] > 1)
                    {
                        continue;
                    }

                    bool grownUpPt(boundaryPoints[pointi] == 1 ? true : false);

                    label nBoundaryNbr = 0;
                    label nGrownUpNbr = 0;
                    label nZoneNbr = 0;
                    forAll(pEdges, pEI)
                    {
                        label edgei = pEdges[pEI];
                        if (cEdgesSet.found(edgei))
                        {
                            edge e = mesh.edges()[edgei];
                            label otherPt(e[0] == pointi ? e[1] : e[0]);

                            if (boundaryPoints[otherPt] == 1)
                            {
                                nGrownUpNbr++;
                            }
                            else if
                            (
                                boundaryPoints[otherPt] != -1
                            )
                            {
                                nBoundaryNbr++;
                            }
                            else if
                            (
                                stationaryPts[otherPt] == -1
                                || stationaryPts[otherPt] == 1
                            )
                            {
                                nZoneNbr++;
                            }

                            if (pointLevel[otherPt] > cLevel)
                            {
                                if (grownUpPt)
                                {
                                    if (boundaryPoints[otherPt] == 1)
                                    {
                                        hPts.append(otherPt);
                                    }
                                }
                                else
                                {
                                    hPts.append(otherPt);
                                }
                            }
                            else
                            {
                                if (grownUpPt)
                                {
                                    if (boundaryPoints[otherPt] == 1)
                                    {
                                        oPts.append(otherPt);
                                    }
                                }
                                else
                                {
                                    oPts.append(otherPt);
                                }
                            }
                        }
                    }

                    if
                    (
                        boundaryPoints[pointi] == -1
                        &&
                        (
                            (
                                nGrownUpNbr > 0
                                && (nBoundaryNbr+nGrownUpNbr) > 1
                            )
                            || nZoneNbr > 1
                        )

                    )
                    {
                        continue;
                    }

                    if (hPts.size() == 4)
                    {
                        label nLayerPts = 0;
                        forAll(hPts, hPI)
                        {
                            label npointi = hPts[hPI];
                            if (layerPointType[npointi] >= 2)
                            {
                                nLayerPts++;
                            }
                        }
                        if (nLayerPts == 0)
                        {
                            hangingNodes[pointi] = true;
                            hangingNbrs[pointi] = hPts;
                        }
                    }
                    else if (oPts.size() == 2)
                    {
                        label nLayerPts = 0;
                        forAll(oPts, oPI)
                        {
                            label npointi = oPts[oPI];
                            if (layerPointType[npointi] >= 2)
                            {
                                nLayerPts++;
                            }
                        }
                        if (nLayerPts == 0)
                        {
                            hangingNodes[pointi] = true;
                            hangingNbrs[pointi] = oPts;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            hangingNodes,
            orEqOp<bool>(),
            false
        );
    }

    scalarField initExpansion(meshPoints.size(), 20);
    pointField outerPos = newPoints;

    for (label iter = 0; iter < nSmoothIter; iter++)
    {
        if (iter % 10 == 0)
        {
            Info<<"Iter : "<<iter<<endl;
        }

        resetPts = -1;

        //Used for carying out stack points and old points field
        // for error reversal
        tPts = newPoints;
        outerPos = newPoints;

        if (debug)
        {
            Time& runTime = const_cast<Time&>(mesh.time());
            runTime++;
            pointField origPts(mesh.points());
            mesh.movePoints(newPoints);
            mesh.write();
            mesh.movePoints(origPts);
        }

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        if (layerParams_.squishTol() > -1 && iter > 10 && (iter % 2 == 0))
        {
            deactivateSquishCells
            (
                layerPointType,
                boundaryPoints,
                newFaceCentres,
                newFaceAreas,
                newCellCentres,
                stationaryPts
            );
        }

        newPoints = vector::zero;

        projectOuterLayerFaces
        (
            iter,
            outerShellPtr,
            layerPointType,
            layerFaceType,
            ftrPointOrigin,
            boundaryPoints,
            stationaryPts,
            stackEdges,
            isMasterFace,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            pointOrigin,
            pointSurfNormal,
            tPts,
            newPoints,
            resetPts
        );

        limitLayerHeight
        (
            ftrPointOrigin,
            layerPointType,
            maxLayerThickness,
            pointSurfNormal,
            pointOrigin,
            newPoints
        );

        resetPoints
        (
            iter,
            controller,
            layerPointType,
            layerFaceType,
            layerEdges,
            medialDist,
            medialLayerHeight,
            pointOrigin,
            pointSurfNormal,
            tPts,
            outerPos,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            newPoints,
            resetPts,
            false
        );

        tPts = newPoints;

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        volumeSmoothPts
        (
            controller,
            isMasterEdge,
            hangingNodes,
            hangingNbrs,
            layerPointType,
            ftrPointOrigin,
            boundaryPoints,
            stationaryPts,
            nLayerCells,
            nNonLayerCells,
            stackEdges,
            featureEdges,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            pointOrigin,
            pointSurfNormal,
            tPts,
            newPoints,
            resetPts
        );

        resetPoints
        (
            iter,
            controller,
            layerPointType,
            layerFaceType,
            layerEdges,
            medialDist,
            medialLayerHeight,
            pointOrigin,
            pointSurfNormal,
            tPts,
            outerPos,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            newPoints,
            resetPts,
            (iter > nSmoothIter/2) // when to correct errors
        );

        tPts = vector::zero;

        reSnap
        (
            grownUpSnapPP,
            grownUpFeatureMeshPtr,
            fzonePP,
            grownUpGeometryPtr,
            grownUpZoneGeometryPtr,
            boundaryPoints,
            stationaryPts,
            newFaceAreas,
            newFaceCentres,
            tPts,
            newPoints
        );

        outerPos = vector::zero;

        reProjectOuter
        (
            ftrPointOrigin,
            layerPointType,
            maxLayerThickness,
            boundaryEdges,
            stackEdges,
            pointSurfNormal,
            pointOrigin,
            newPoints,
            outerPos
        );

        stretchLayerMesh
        (
            controller,
            layerPointType,
            boundaryEdges,
            stackEdges,
            boundaryPoints,
            stationaryPts,
            pp,
            meshPoints,
            layerMethod,
            totalNumLayers,
            layerCount,
            pointSurfNormal,
            pointOrigin,
            unsmoothedPointNormals,
            outerPos,
            initExpansion,
            newPoints,
            tPts
        );

        if (iter == nSmoothIter-1)
        {
            if (controller.algorithm() == meshControl::EXTRUDE)
            {
                correctConcaveCells
                (
                    layerPointType,
                    layerFaceType,
                    newFaceCentres,
                    newCellCentres,
                    newFaceAreas,
                    newCellVolumes,
                    tPts,
                    newPoints
                );
            }

            if (adjustFinalNLayers)
            {
                adjustNumLayers
                (
                    pp,
                    layerMethod,
                    initExpansion,
                    pointOrigin,
                    outerPos,
                    totalNumLayers,
                    adjustedNLayers
                );
            }
        }

        if (debug && outerShellPtr.valid())
        {
            simpleVTKWriter
            (
                outerShellPtr()(),
                newPoints
            ).write("outerShell"+Foam::name(iter)+".vtk");
        }
    }

    //reuse tPts field for calculating point displacement
    tPts = (newPoints-mesh.points());

    syncTools::syncPointList
    (
        mesh,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = mesh.points() + tPts;

    mesh.movePoints(newPoints);
}


void Foam::layerManipulate::stretchLayers
(
    const meshControl& controller
)
{
    Info<<"Stretching layer cells "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const labelIOList& layerCells = meshRefiner_.layerCells()();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    boolList layerEdges(mesh.nEdges(), false);
    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            const labelList& cellEdges = mesh.cellEdges()[celli];
            forAll(cellEdges, cEI)
            {
                layerEdges[cellEdges[cEI]] = true;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    //Set layer growth method
    labelList ppPtLevel(meshPoints.size(), -1);
    labelList ppPtLevelMin(meshPoints.size(), labelMax);
    forAll(meshPoints, pti)
    {
        const labelList& pFaces = pp.pointFaces()[pti];
        forAll(pFaces, pfi)
        {
            label meshFaceI = pp.addressing()[pFaces[pfi]];
            label own = mesh.faceOwner()[meshFaceI];
            ppPtLevel[pti] = max(cellLevel[own],ppPtLevel[pti]);
            ppPtLevelMin[pti] = min(cellLevel[own],ppPtLevelMin[pti]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        ppPtLevel,
        maxEqOp<label>(),
        label(-1)         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        ppPtLevelMin,
        minEqOp<label>(),
        labelMax         // null value
    );

    List<Tuple2<word,scalar>> layerMethod(meshPoints.size());
    scalarField ppMaxLayerThickness(meshPoints.size(), GREAT);
    setLayerMethod
    (
        pp,
        meshPoints,
        ppPtLevel,
        ppPtLevelMin,
        layerMethod,
        ppMaxLayerThickness
    );

    boolList boundaryFaces(mesh.nFaces(), false);
    boolList boundaryEdges(mesh.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    autoPtr<indirectPrimitivePatch> grownUpPPPtr = makeGrownUpPatch();
    indirectPrimitivePatch& grownUpPP = grownUpPPPtr();

    autoPtr<indirectPrimitivePatch> fzonePPPtr = makeFaceZonePatch();
    indirectPrimitivePatch& fzonePP = fzonePPPtr();

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
    boolList featurePts(mesh.nPoints(), false);
    boolList featureEdges(mesh.nEdges(), false);

    //Mark zone feature edges
    if
    (
        layerParams_.dualReSnapZones()
        && returnReduce(fzonePP.size(), sumOp<label>()) != 0
    )
    {
        const labelList zoneEdges
        (
            fzonePP.meshEdges
            (
                mesh.edges(),
                mesh.pointEdges()
            )
        );

        autoPtr<boolList> flipMap
        (
            new boolList(fzonePP.size(), false)
        );
        labelList neiCellZone;
        labelList cellToZone(mesh.nCells(), -1);
        forAll(mesh.cells(), celli)
        {
            label zonei = mesh.cellZones().whichZone(celli);
            cellToZone[celli] = zonei;
        }
        syncTools::swapBoundaryCellList(mesh, cellToZone, neiCellZone);

        forAll(fzonePP, i)
        {
            label facei = fzonePP.addressing()[i];
            label ownZoneI = mesh.cellZones().whichZone
            (
                mesh.faceOwner()[facei]
            );

            label  neiZone = -1;
            if (mesh.isInternalFace(facei))
            {
                neiZone = mesh.cellZones().whichZone
                (
                    mesh.faceNeighbour()[facei]
                );
            }
            else
            {
                neiZone = neiCellZone[facei-mesh.nInternalFaces()];
            }
            if (ownZoneI < neiZone)
            {
                flipMap()[i] = true;
            }
        }

        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            fzonePP,
            zoneEdges,
            0.8191,
            0.8191,
            -GREAT,
            flipMap
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(zoneEdges, edgei)
        {
            if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
                || eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = zoneEdges[edgei];
                edge e = mesh.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }
    }

    //Calculate feature edges grown up patch
    {
        const labelList meshEdgesGrownUp
        (
            grownUpPP.meshEdges
            (
                mesh.edges(),
                mesh.pointEdges()
            )
        );

        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            grownUpPP,
            meshEdgesGrownUp,
            0.8191,
            0.8191
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(meshEdgesGrownUp, edgei)
        {
            if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
                || eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = meshEdgesGrownUp[edgei];
                edge e = mesh.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            featureEdges,
            orEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh,
            featurePts,
            orEqOp<bool>(),
            false
        );

        labelList nFeatureEdges(mesh.nPoints(), 0);
        forAll(mesh.points(), pointi)
        {
            labelList pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                if (isMasterEdge[pEdges[pEI]] && featureEdges[pEdges[pEI]])
                {
                    nFeatureEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nFeatureEdges,
            plusEqOp<label>(),
            label(0)
        );

        forAll(mesh.points(), pointi)
        {
            if (boundaryPoints[pointi] == -1)
            {
                if (featurePts[pointi])
                {
                    if (nFeatureEdges[pointi] > 2)
                    {
                        boundaryPoints[pointi] = 3;
                    }
                    else
                    {
                        boundaryPoints[pointi] = 2;
                    }
                }
            }
        }

        forAll(grownUpPP.meshPoints(), ptI)
        {
            label pointi = grownUpPP.meshPoints()[ptI];

            if (boundaryPoints[pointi] == -1)
            {
                boundaryPoints[pointi] = 1;
            }
        }
    }

    forAll(pp.edges(), edgei)
    {
        boundaryEdges[meshEdges[edgei]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh.nCells(), 0);
    labelList stackEdgeLayer(mesh.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh.nFaces(), -1);
    labelList layerPointType(mesh.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    boolList stackEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        if (stackEdgeLayer[edgei] != 0)
        {
            stackEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    pointField pointNormals,unsmoothedPointNormals;
    //Calculate area weighted point normals (feature points fixed)
    {
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            0.707,//-0.5,
            0.707//-0.5
        );

        boolList excludedFaces(pp.size(), false);
        pointNormals = eClass.calculatePointNormals
        (
            excludedFaces,
            layerParams_.nSmoothNormalsExtrude(),
            true
        );
        unsmoothedPointNormals =
            eClass.calculatePointNormals(excludedFaces, 0, true);
    }

    mesh.clearOut();

    labelList layerCount(mesh.nPoints(), -1);

    pointField pointOrigin(mesh.nPoints(), vector::zero);
    pointField pointSurfNormal(mesh.nPoints(), vector::zero);

    scalarField maxLayerThickness(mesh.nPoints(), GREAT);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = label(0);
        pointOrigin[meshPointI] = pp.localPoints()[ptI];
        pointSurfNormal[meshPointI] = -pointNormals[ptI];
        maxLayerThickness[meshPointI] = ppMaxLayerThickness[ptI];
    }

    labelList nVisited(mesh.nPoints(), 0);

    boolList edgeSet(mesh.nEdges(), false);

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            pointOrigin,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointSurfNormal,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncEdgeList
        (
            mesh,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        forAll(mesh.edges(), edgei)
        {
            if (!isMasterEdge[edgei])
            {
                continue;
            }

            if (!boundaryEdges[edgei] && stackEdges[edgei] && !edgeSet[edgei])
            {
                edge e = mesh.edges()[edgei];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgei] = true;

                        nVisited[e[1]]++;

                        maxLayerThickness[e[1]] = maxLayerThickness[e[0]];
                        pointOrigin[e[1]] += pointOrigin[e[0]];
                        pointSurfNormal[e[1]] += pointSurfNormal[e[0]];
                        nSet++;
                    }
                }

                if (!edgeSet[edgei] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgei] = true;

                        nVisited[e[0]]++;

                        maxLayerThickness[e[0]] = maxLayerThickness[e[1]];
                        pointOrigin[e[0]] += pointOrigin[e[1]];
                        pointSurfNormal[e[0]] += pointSurfNormal[e[1]];
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nVisited,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.points(), pointi)
    {
        if (nVisited[pointi] > 1)
        {
            pointOrigin[pointi] /= nVisited[pointi];
            pointSurfNormal[pointi] /= nVisited[pointi];
        }
    }

    //Check for unitinitilised hanging points
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh.points(), pointi)
        {
            if (layerPointType[pointi] == 2 && nVisited[pointi] == 0)
            {
                const labelList& pEdges = mesh.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    edge e = mesh.edges()[edgei];
                    label otherPt(e[0] ==  pointi ? e[1] : e[0]);
                    if (layerEdges[edgei] && layerPointType[otherPt] == 2)
                    {
                        maxLayerThickness[pointi] = min
                        (
                            maxLayerThickness[otherPt],
                            maxLayerThickness[pointi]
                        );
                        pointOrigin[pointi] += pointOrigin[otherPt];
                        pointSurfNormal[pointi] += pointSurfNormal[otherPt];
                        nVisited[pointi]++;
                    }
                }
            }
            else
            {
                nVisited[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nVisited,
            plusEqOp<label>(),
            label(0)
        );

        syncTools::syncPointList
        (
            mesh,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        syncTools::syncPointList
        (
            mesh,
            pointOrigin,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointSurfNormal,
            plusEqOp<vector>(),
            vector::zero
        );

        forAll(mesh.points(), pointi)
        {
            if (nVisited[pointi] > 1)
            {
                pointOrigin[pointi] /= nVisited[pointi];
                pointSurfNormal[pointi] /= nVisited[pointi];
            }
        }
    }

    labelList totalNumLayers(mesh.nPoints(), -1);
    forAll(mesh.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            face f = mesh.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                totalNumLayers[pointi] = max
                (
                    totalNumLayers[pointi],layerCount[pointi]
                );
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(mesh.edges(), edgei)
        {
            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }


    //Prevent faceZones from having blended profile
    labelList stationaryPts(mesh.nPoints(), -2);
    forAll(fzonePP.meshPoints(), pointi)
    {
        label meshPointI = fzonePP.meshPoints()[pointi];
        stationaryPts[meshPointI] = -1;
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    scalarField initExpansion(meshPoints.size(), 20);
    pointField newPoints(mesh.points());
    pointField tPts(mesh.points());

    pointField outerPos(mesh.points());
    boolList pointVisited(mesh.nPoints(), false);
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            pointVisited[pointi] = true;
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            outerPos,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointVisited,
            orEqOp<bool>(),
            false
        );

        forAll(mesh.edges(), edgei)
        {
            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                if
                (
                    pointVisited[e[0]] && !pointVisited[e[1]]
                )
                {
                    pointVisited[e[1]] = true;
                    outerPos[e[1]] = outerPos[e[0]];
                    nSet++;
                }
                else if
                (
                    pointVisited[e[1]] && !pointVisited[e[0]]
                )
                {
                    pointVisited[e[0]] = true;
                    outerPos[e[0]] = outerPos[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    stretchLayerMesh
    (
        controller,
        layerPointType,
        boundaryEdges,
        stackEdges,
        boundaryPoints,
        stationaryPts,
        pp,
        meshPoints,
        layerMethod,
        totalNumLayers,
        layerCount,
        pointSurfNormal,
        pointOrigin,
        unsmoothedPointNormals,
        outerPos,
        initExpansion,
        newPoints,
        tPts
    );

    mesh.movePoints(newPoints);
}


void Foam::layerManipulate::volumeSmoothPts
(
    const meshControl& controller,
    const PackedBoolList& isMasterEdge,
    const boolList& hangingNodes,
    const List<labelList>& hangingNbrs,
    const labelList& layerPointType,
    const labelList& ftrPointOrigin,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const scalarField& nLayerCells,
    const scalarField& nNonLayerCells,
    const boolList& stackEdges,
    const boolList& featureEdges,
    const scalarField& newCellVolumes,
    const vectorField& newFaceAreas,
    const pointField& newFaceCentres,
    const pointField& newCellCentres,
    const pointField& pointOrigin,
    const vectorField& pointSurfNormal,
    const pointField& tPts,
    pointField& newPoints,
    labelList& resetPts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const scalar edge0Length = meshRefiner_.meshCutter().level0EdgeLength();

    scalar weight = 0.5;
    scalar invWeight = 1-weight;
    scalar reducedWeight = 0.1;
    scalar invReducedWeight = 1-reducedWeight;

    scalarField sumVol(mesh.nPoints(), 0);
    pointField avePt(mesh.nPoints(), 0);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelIOList& layerCells = meshRefiner_.layerCells()();

    scalarField nonLayerVolWeights(mesh.nPoints(),0);
    scalarField layerVolWeights(mesh.nPoints(),0);
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            if (boundaryPoints[pointi] == 1)
            {
                const labelList& pFaces = mesh.pointFaces()[pointi];
                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];
                    label patchi = patches.whichPatch(facei);
                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        label celli = mesh.faceOwner()[facei];
                        label  level = cellLevel[celli];
                        scalar len = edge0Length / pow(2., level);
                        scalar faceWeight = sqrt(mag(newFaceAreas[facei]))/len;
                        if (layerCells[celli] > -1)
                        {
                            layerVolWeights[pointi] += faceWeight;
                        }
                        else
                        {
                            nonLayerVolWeights[pointi] += faceWeight;
                        }
                    }
                }
            }
            else
            {
                const labelList& pCells = mesh.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    label  level = cellLevel[celli];
                    scalar len = edge0Length / pow(2., level);
                    scalar volWeight = cbrt(mag(newCellVolumes[celli]))/len;

                    if (layerCells[celli] > -1)
                    {
                        layerVolWeights[pointi] += volWeight;
                    }
                    else
                    {
                        nonLayerVolWeights[pointi] += volWeight;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        layerVolWeights,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nonLayerVolWeights,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    //Calculate weights at interface of layer and non-layer cells
    scalar nliw = layerParams_.layerInterfaceWeights();
    scalar liw = scalar(1) - nliw;
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            scalar nCells =  nNonLayerCells[pointi]
                + nLayerCells[pointi];
            scalar nVolWeights = nonLayerVolWeights[pointi]
                + layerVolWeights[pointi];
            if (nVolWeights > 0)
            {
                scalar aveNLW = nonLayerVolWeights[pointi]
                    / nNonLayerCells[pointi];
                scalar aveLW = layerVolWeights[pointi]
                    / nLayerCells[pointi];
                nonLayerVolWeights[pointi] =
                   nliw*(aveNLW/nNonLayerCells[pointi]);
                layerVolWeights[pointi] = liw*(aveLW/nLayerCells[pointi]);
            }
            else if (nCells > 0)
            {
                nonLayerVolWeights[pointi] = 1.0 / nCells;
                layerVolWeights[pointi] = 1.0 / nCells;
            }
        }
    }

    forAll(newPoints, pointi)
    {
        if
        (
            boundaryPoints[pointi] != 0
            && (stationaryPts[pointi] < 0 || stationaryPts[pointi] == 3)
            && (layerPointType[pointi] == -1 || layerPointType[pointi] >= 2)
         )
        {
            if (boundaryPoints[pointi] > 0)
            {
                if
                (
                    layerPointType[pointi] == 2
                    && boundaryPoints[pointi] == 1
                )
                {
                    labelList pFaces = mesh.pointFaces()[pointi];
                    forAll(pFaces, pFI)
                    {
                        label facei = pFaces[pFI];
                        label patchI = patches.whichPatch(facei);

                        if (patchI != -1 && !patches[patchI].coupled())
                        {
                            label own = mesh.faceOwner()[facei];
                            if (layerCells[own] > -1)
                            {
                                scalar wf = -1;
                                if
                                (
                                    controller.algorithm()
                                    == meshControl::DUAL
                                 )
                                {
                                    wf = layerVolWeights[pointi]
                                        *mag(newFaceAreas[facei]);
                                }
                                else
                                {
                                    label  level =
                                        cellLevel[mesh.faceOwner()[facei]];
                                    scalar len =
                                        edge0Length / pow(2., level);
                                    wf = layerVolWeights[pointi]
                                        *sqrt(mag(newFaceAreas[facei]))/len;
                                }
                                avePt[pointi] +=
                                    (wf * newFaceCentres[facei]);
                                sumVol[pointi] += wf;
                            }
                            else
                            {
                                scalar wf = -1;
                                if
                                (
                                    controller.algorithm()
                                    == meshControl::DUAL
                                 )
                                {
                                    wf = nonLayerVolWeights[pointi]
                                        *mag(newFaceAreas[facei]);
                                }
                                else
                                {
                                    label level =
                                        cellLevel[mesh.faceOwner()[facei]];
                                    scalar len = edge0Length
                                        / pow(2., level);
                                    wf = nonLayerVolWeights[pointi]
                                        *sqrt(mag(newFaceAreas[facei]))/len;
                                }
                                avePt[pointi] +=
                                    (wf * newFaceCentres[facei]);
                                sumVol[pointi] += wf;
                            }
                        }
                    }
                }
                else
                {
                    if (boundaryPoints[pointi] == 1)
                    {
                        if (hangingNodes[pointi])
                        {
                            const labelList& nbrs = hangingNbrs[pointi];
                            forAll(nbrs, nbri)
                            {
                                avePt[pointi] += newPoints[nbrs[nbri]];
                                sumVol[pointi] += scalar(1);
                            }
                        }
                        else
                        {
                            labelList pFaces = mesh.pointFaces()[pointi];
                            forAll(pFaces, pFI)
                            {
                                label facei = pFaces[pFI];
                                label patchI = patches.whichPatch(facei);
                                if
                                (
                                    patchI != -1
                                    && !patches[patchI].coupled()
                                 )
                                {
                                    scalar wf = -1;
                                    if
                                    (
                                        controller.algorithm()
                                        == meshControl::DUAL
                                     )
                                    {
                                        wf = mag(newFaceAreas[facei]);
                                    }
                                    else
                                    {
                                        label  level = cellLevel
                                            [mesh.faceOwner()[facei]];
                                        scalar len = edge0Length
                                            / pow(2., level);
                                        wf = sqrt(mag(newFaceAreas[facei]))
                                            /len;
                                    }

                                    avePt[pointi] +=
                                        (wf * newFaceCentres[facei]);
                                    sumVol[pointi] += wf;
                                }
                            }
                        }
                    }
                    else if (boundaryPoints[pointi] == 2)
                    {
                        labelList pEdges = mesh.pointEdges()[pointi];
                        forAll(pEdges, pEI)
                        {
                            label edgei = pEdges[pEI];
                            if (featureEdges[edgei] && isMasterEdge[edgei])
                            {
                                edge e = mesh.edges()[edgei];
                                scalar edgeWeight = e.mag(newPoints);
                                if
                                (
                                    controller.algorithm()
                                    != meshControl::DUAL
                                 )
                                {
                                    label level = max
                                    (
                                        pointLevel[e[0]],
                                        pointLevel[e[1]]
                                    );
                                    scalar len = edge0Length
                                        / pow(2., level);
                                    edgeWeight /= len;
                                }
                                avePt[pointi] +=
                                    (edgeWeight * e.centre(newPoints));
                                sumVol[pointi] += edgeWeight;
                            }
                        }
                    }
                }
            }
            else
            {
                if (hangingNodes[pointi] && layerPointType[pointi] != 2)
                {
                    const labelList& nbrs = hangingNbrs[pointi];
                    forAll(nbrs, nbri)
                    {
                        avePt[pointi] += newPoints[nbrs[nbri]];
                        sumVol[pointi] += scalar(1);
                    }
                }
                else
                {
                    labelList pCells = mesh.pointCells()[pointi];
                    forAll(pCells, pCI)
                    {
                        label celli = pCells[pCI];
                        if (layerPointType[pointi] == 2)
                        {
                            scalar volWeight = -1;
                            if (controller.algorithm() == meshControl::DUAL)
                            {
                                volWeight = mag(newCellVolumes[celli]);
                            }
                            else
                            {
                                label  level = cellLevel[celli];
                                scalar len = edge0Length / pow(2., level);
                                volWeight =
                                    cbrt(mag(newCellVolumes[celli]))/len;
                            }

                            if (layerCells[celli] > -1)
                            {
                                scalar wf = layerVolWeights[pointi]
                                    * volWeight;
                                avePt[pointi] += (wf*newCellCentres[celli]);
                                sumVol[pointi] += wf;
                            }
                            else
                            {
                                scalar wf = nonLayerVolWeights[pointi]
                                    * volWeight;
                                avePt[pointi] += (wf*newCellCentres[celli]);
                                sumVol[pointi] += wf;
                            }
                        }
                        else
                        {
                            scalar volWeight = -1;
                            if (controller.algorithm() == meshControl::DUAL)
                            {
                                volWeight = mag(newCellVolumes[celli]);
                            }
                            else
                            {
                                label  level = cellLevel[celli];
                                scalar len = edge0Length / pow(2., level);
                                volWeight =
                                    cbrt(mag(newCellVolumes[celli]))/len;
                            }
                            avePt[pointi] +=
                                (volWeight * newCellCentres[celli]);
                            sumVol[pointi] += volWeight;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        sumVol,
        plusEqOp<scalar>(),
        scalar(0)
     );

    syncTools::syncPointList
    (
        mesh,
        avePt,
        plusEqOp<point>(),
        vector::zero
     );

    resetPts = -1;

    forAll(newPoints, pointi)
    {
        if
        (
            boundaryPoints[pointi] != 0
            && (layerPointType[pointi] == -1 || layerPointType[pointi] >= 2)
        )
        {
            if (sumVol[pointi] > VSMALL)
            {
                if (stationaryPts[pointi] < 0 && layerPointType[pointi] >= 2)
                {
                    point updatedPt;

                    if (ftrPointOrigin[pointi] > 0)
                    {
                        updatedPt = invReducedWeight*tPts[pointi]
                            +reducedWeight*(avePt[pointi]/sumVol[pointi]);
                    }
                    else
                    {
                        updatedPt = invWeight*tPts[pointi]
                            +weight*(avePt[pointi]/sumVol[pointi]);
                    }

                    newPoints[pointi] = updatedPt;

                    //Do not freeze boundary points
                    if (boundaryPoints[pointi] > 0)
                    {
                        continue;
                    }

                    labelList pEdges = mesh.pointEdges()[pointi];

                    forAll(pEdges, pEI)
                    {
                        label edgei = pEdges[pEI];

                        if (stackEdges[edgei])
                        {
                            edge e = mesh.edges()[edgei];
                            label otherPt = (e[0] ==  pointi ? e[1] : e[0]);

                            vector toSurf = updatedPt
                                - pointOrigin[otherPt];
                            toSurf /= (mag(toSurf) + VSMALL);
                            vector sNorm = pointSurfNormal[otherPt];
                            sNorm /= (mag(sNorm) + VSMALL);

                            vector toSurfOrig =
                                tPts[pointi] - pointOrigin[otherPt];
                            toSurfOrig /= (mag(toSurfOrig) + VSMALL);

                            scalar dP1 = (toSurf&sNorm);
                            scalar dP2 = (toSurfOrig&sNorm);

                            if
                            (
                                dP1 < dP2 && dP1 < 0.866
                                && ftrPointOrigin[otherPt] != -1
                             )
                            {
                                resetPts[pointi] = label(1);
                                break;
                            }
                        }
                    }
                }
                else
                {
                    if (stationaryPts[pointi] > -1)
                    {
                        if (stationaryPts[pointi] == 3)
                        {
                            newPoints[pointi] = invReducedWeight*tPts[pointi]
                                + reducedWeight*(avePt[pointi]/sumVol[pointi]);
                        }
                        else
                        {
                            newPoints[pointi] = tPts[pointi];
                        }
                    }
                    else
                    {
                        newPoints[pointi] = invWeight*tPts[pointi]
                            +weight*(avePt[pointi]/sumVol[pointi]);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        resetPts,
        maxEqOp<label>(),
        label(-1)
     );
}


void Foam::layerManipulate::stretchLayerMesh
(
    const meshControl& controller,
    const labelList& layerPointType,
    const boolList& boundaryEdges,
    const boolList& stackEdges,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const indirectPrimitivePatch& pp,
    const labelList& meshPoints,
    const List<Tuple2<word,scalar>>& layerMethod,
    const labelList& totalNumLayers,
    const labelList& layerCount,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    const pointField& unsmoothedPointNormals,
    const pointField& outerPos,
    scalarField& initExpansion,
    pointField& newPoints,
    pointField& tPts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    scalar tol = 0.0001;
    scalar maxExpansion = 4;
    scalar limitExpansion = 20;
    scalar minStretch = layerParams_.minStretch();
    point greatPoint(GREAT, GREAT, GREAT);

    scalarField expansionRatio(mesh.nPoints(), -1);

    const List<wordList>& layerSpec = layerParams_.layerSpec();

    const label fixedLayerCount = layerParams_.nFCHLayers();
    const List<Switch> fixedFCH = layerParams_.fixedFCH();
    bool setExactFCH = false;
    forAll(fixedFCH, patchi)
    {
        if (fixedFCH[patchi])
        {
            word method = layerSpec[patchi][1];
            if (method == "fch" || method == "rfch")
            {
                setExactFCH = true;
                break;
            }
        }
    }

    vectorField exactFCH;
    boolList fixedFCHPts;
    if (setExactFCH)
    {
        exactFCH.setSize(mesh.nPoints(), vector::zero);
        fixedFCHPts.setSize(meshPoints.size(), false);
        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label patchi = patches.whichPatch(meshFaceI);
            if (fixedFCH[patchi])
            {
                face f = pp.localFaces()[i];
                forAll(f,fp)
                {
                    fixedFCHPts[f[fp]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            fixedFCHPts,
            orEqOp<bool>(),
            false         // null value
        );
    }

    // Allow monotonic smoothing of fch field to handle patch discontinuities
    autoPtr<scalarField> smoothedFCHPtr;
    bool fchSmoothActive = false;
    if (layerParams_.nFCHSmooth() > 0)
    {
        fchSmoothActive = true;
        smoothedFCHPtr = new scalarField(meshPoints.size(),scalar(-1));
        scalarField& smoothedFCH = smoothedFCHPtr();
        forAll(meshPoints, pointi)
        {
            if (layerMethod[pointi].first() == "fch")
            {
                label meshpointi = meshPoints[pointi];
                scalar firstCellHeight = layerMethod[pointi].second();
                vector layerVec = outerPos[meshpointi]-pointOrigin[meshpointi];
                scalar layerHeight = mag(layerVec);
                if (layerHeight > VSMALL)
                {
                    vector unitLayerVec = layerVec / layerHeight;
                    scalar dProd = max
                    (
                        0.5,
                        mag(unitLayerVec & unsmoothedPointNormals[pointi])
                    );
                    //correct first cell height
                    firstCellHeight /= dProd;
                }
                smoothedFCH[pointi] = firstCellHeight;
            }
        }

        boolList smoothedFaces(pp.size(), false);
        forAll(pp,i)
        {
            const face& lf = pp.localFaces()[i];
            label nfch = 0;
            forAll(lf,lfp)
            {
                if (layerMethod[lf[lfp]].first() == "fch")
                {
                    nfch++;
                }
            }
            if (nfch == lf.size())
            {
                smoothedFaces[i] = true;
            }
        }

        boolList smoothedPts(meshPoints.size(), false);
        labelList sumNFaces(meshPoints.size(), label(0));
        forAll(meshPoints, pointi)
        {
            label meshpointi = meshPoints[pointi];
            if (boundaryPoints[meshpointi] != 1)
            {
                const labelList& pFaces = pp.pointFaces()[pointi];
                label nSmoothedFaces = 0;
                forAll(pFaces, pfi)
                {
                    label facei = pFaces[pfi];
                    if (smoothedFaces[facei])
                    {
                        sumNFaces[pointi]++;
                        nSmoothedFaces++;
                    }
                }
                if (nSmoothedFaces == pFaces.size())
                {
                    smoothedPts[pointi] = true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            smoothedPts,
            orEqOp<bool>(),
            false         // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            sumNFaces,
            plusEqOp<label>(),
            label(0)        // null value
        );

        scalarField sumFCH(meshPoints.size(), 0.0);
        for
        (
            label nSmooth = 0; nSmooth < layerParams_.nFCHSmooth(); nSmooth++
        )
        {
            sumFCH = 0;
            forAll(meshPoints, pointi)
            {
                if (smoothedPts[pointi])
                {
                    const labelList& pFaces = pp.pointFaces()[pointi];
                    forAll(pFaces, pfi)
                    {
                        label facei = pFaces[pfi];
                        scalar aveFCH = scalar(0);
                        if (smoothedFaces[facei])
                        {
                            const face& lf = pp.localFaces()[facei];
                            forAll(lf,lpi)
                            {
                                label lpointi = lf[lpi];
                                scalar sfch = smoothedFCH[lpointi];
                                aveFCH += sfch;
                            }
                            aveFCH /= lf.size();
                            sumFCH[pointi] += aveFCH;
                       }
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                meshPoints,
                sumFCH,
                plusEqOp<scalar>(),
                scalar(0)        // null value
            );

            forAll(meshPoints, pointi)
            {
                if (smoothedPts[pointi] && sumNFaces[pointi] > 0)
                {
                    scalar sFCH = sumFCH[pointi] / sumNFaces[pointi];
                    if (sFCH < smoothedFCH[pointi])
                    {
                        smoothedFCH[pointi] = 0.5*(smoothedFCH[pointi] + sFCH);
                    }
                }
            }
        }
    }

    //Calculate stretching at primitive patch points
    forAll(meshPoints, ptI)
    {
        label meshPointI = meshPoints[ptI];
        if (layerMethod[ptI].first() == "fch")
        {
            label numLayers = totalNumLayers[meshPointI];
            vector layerVec = outerPos[meshPointI]-pointOrigin[meshPointI];
            scalar layerHeight = mag(layerVec);

            if (layerHeight < VSMALL)
            {
                expansionRatio[meshPointI] = 1.0;
                continue;
            }

            scalar firstCellHeight = layerMethod[ptI].second();
            if (fchSmoothActive)
            {
                firstCellHeight = smoothedFCHPtr()[ptI];
            }
            else
            {
                vector unitLayerVec = layerVec / layerHeight;
                scalar dProd = max
                (
                    0.5,
                    mag(unitLayerVec & unsmoothedPointNormals[ptI])
                );
                //correct first cell height
                firstCellHeight /= dProd;
            }

            if (firstCellHeight >= layerHeight)
            {
                expansionRatio[meshPointI] = 1.0;
                continue;
            }

            scalar xn = initExpansion[ptI];
            scalar x = -1.0;
            label repeat = 0;
            if (xn < 1.0 + tol && xn > 1.0 - tol)
            {
                xn = maxExpansion;
                x = -1.0;
                repeat = 0;
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;
                    if
                    (
                        f2 < VSMALL && f2 > -VSMALL
                    )
                    {
                        WarningInFunction
                            << "f2 is less than SMALL : " << f2
                            << " fch : " << firstCellHeight
                            << " number layers : " << numLayers
                            << " layerHeight : " << layerHeight
                            << " expansion : " << x
                            << " at location : " << pointOrigin[meshPointI]
                            << endl;
                        xn =1.0;
                        break;
                    }
                    xn = x - f1/f2;
                    xn  = max(scalar(0), min(xn,limitExpansion));
                }
            }
            else
            {
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;

                    if
                    (
                        f2 < VSMALL && f2 > -VSMALL
                    )
                    {
                        WarningInFunction
                            << "f2 is less than SMALL : " << f2
                            << " fch : " << firstCellHeight
                            << " number layers : " << numLayers
                            << " layerHeight : " << layerHeight
                            << " expansion : " << x
                            << " at location : " << pointOrigin[meshPointI]
                            << endl;
                        xn =1.0;
                        break;
                    }
                    xn = x - f1/f2;
                    xn  = max(scalar(0), min(xn,limitExpansion));
                }
            }

            if (xn < 1.0 + tol && xn > 1.0 - tol)
            {
                xn = 0.0;
                x = -1.0;
                repeat = 0;
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;
                    if (f2 < VSMALL && f2 > -VSMALL)
                    {
                        WarningInFunction
                            << "f2 is less than SMALL : " << f2
                            << " fch : " << firstCellHeight
                            << " number layers : " << numLayers
                            << " layerHeight : " << layerHeight
                            << " expansion : " << x
                            << " at location : " << pointOrigin[meshPointI]
                            << endl;
                        xn =1.0;
                        break;
                    }
                    xn = x - f1/f2;
                    xn  = max(scalar(0), min(xn,limitExpansion));
                }
            }

            scalar limitMin = max(minStretch,xn);
            initExpansion[ptI] = limitMin;
            expansionRatio[meshPointI] = min(limitMin,maxExpansion);
            if (setExactFCH && fixedFCHPts[ptI])
            {
                xn = expansionRatio[meshPointI];
                scalar actualFCH;
                scalar layer2height;
                if (xn < 1.0 - tol || xn > 1.0 + tol)
                {
                    actualFCH = layerHeight*(1. - xn)
                        /(1.-pow(xn, numLayers));
                    layer2height = (1.+xn)*actualFCH;
                }
                else
                {
                    actualFCH = layerHeight/numLayers;
                    layer2height = 2*actualFCH;
                }
                scalar inputFCH = fchSmoothActive ?
                    smoothedFCHPtr()[ptI] : layerMethod[ptI].second();

                if (inputFCH < 0.8*layer2height)
                {
                    exactFCH[meshPointI] =
                        -inputFCH*unsmoothedPointNormals[ptI];
                }
            }
        }
        else
        {
            expansionRatio[meshPointI] = layerMethod[ptI].second();
        }
    }

    tPts = -greatPoint;

    //Check if layer blending required
    bool blend = false;
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        blend = layerParams_.extrudeBlend();
    }

    //Calculate layer stretched layer profile
    {
        boolList pointVisited(mesh.nPoints(), false);
        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            pointVisited[meshPointI] = true;
        }

        forAll(mesh.points(), pointi)
        {
            if (layerPointType[pointi] >= 2)
            {
                tPts[pointi] = outerPos[pointi];
            }
        }

        while (true)
        {
            label nSet = 0;

            syncTools::syncPointList
            (
                mesh,
                expansionRatio,
                maxEqOp<scalar>(),
                maxExpansion
            );

            if (setExactFCH)
            {
                syncTools::syncPointList
                (
                    mesh,
                    exactFCH,
                    maxMagSqrEqOp<point>(),
                    vector::zero
                );
            }

            syncTools::syncPointList
            (
                mesh,
                pointVisited,
                orEqOp<bool>(),
                false
             );

            forAll(mesh.edges(), edgei)
            {
                if (!boundaryEdges[edgei] && stackEdges[edgei])
                {
                    edge e = mesh.edges()[edgei];

                    label pt0 = -1;
                    label pt1 = -1;

                    if (pointVisited[e[0]] && !pointVisited[e[1]])
                    {
                        pt0 = e[0];
                        pt1 = e[1];
                    }
                    else if (pointVisited[e[1]] && !pointVisited[e[0]])
                    {
                        pt0 = e[1];
                        pt1 = e[0];
                    }

                    if
                    (
                        pt0 == -1
                        || boundaryPoints[pt1] == 0
                        || layerPointType[pt1] >=2
                    )
                    {
                        continue;
                    }

                    pointVisited[pt1] = true;
                    expansionRatio[pt1] = expansionRatio[pt0];
                    if (setExactFCH)
                    {
                        exactFCH[pt1] = exactFCH[pt0];
                    }
                    vector dir = (outerPos[pt1]-pointOrigin[pt1]);
                    scalar str = expansionRatio[pt0];

                    label lcount = layerCount[pt1];

                    if
                    (
                        setExactFCH && lcount <= fixedLayerCount
                        && mag(exactFCH[pt0]) > 0
                    )
                    {
                        tPts[pt1] = pointOrigin[pt1] + lcount*exactFCH[pt0];
                    }
                    else if (str > 1.0 + tol || str < 1.0 - tol)
                    {
                        scalar lh = mag(dir) + VSMALL;
                        scalar fch = lh *
                        (
                            (1.-str)/(1.-pow(str,totalNumLayers[pt1]))
                        );
                        dir /= lh;
                        scalar ratio = fch*(1.-pow(str,lcount))
                            /(1.-str);

                        if
                        (
                            blend && boundaryPoints[pt1] != 1
                            && stationaryPts[pt1] != -1
                        )
                        {
                            scalar weight = min
                            (
                                scalar(1.0),
                                max
                                (
                                    scalar(0.0),ratio/lh
                                )
                            );

                            weight = sin
                            (
                                weight
                                * 0.5*Foam::constant::mathematical::pi
                            );
                            weight = Foam::pow(weight, scalar(0.25));

                            dir = (1.-weight)*pointSurfNormal[pt1]
                                + weight*dir;
                        }
                        tPts[pt1] = pointOrigin[pt1]
                            + ratio*dir;
                    }
                    else
                    {
                        scalar ratio = lcount / scalar(totalNumLayers[pt1]);
                        if
                        (
                            blend && boundaryPoints[pt1] != 1
                            && stationaryPts[pt1] != -1
                        )
                        {
                            scalar lh = mag(dir) + VSMALL;
                            vector udir(dir/lh);

                            scalar weight = min
                            (
                                scalar(1.0),
                                max
                                (
                                    scalar(0.0),ratio
                                )
                            );

                            weight = sin
                            (
                                weight
                                * 0.5*Foam::constant::mathematical::pi
                            );
                            weight = Foam::pow(weight, scalar(0.25));

                            dir = lh*
                            (
                                (1.-weight)*pointSurfNormal[pt1]
                                + weight*udir
                            );
                        }
                        tPts[pt1] = pointOrigin[pt1]
                            + ratio*dir;
                    }
                    nSet++;
                }
            }

            if (returnReduce(nSet, sumOp<label>()) == 0)
            {
                break;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        tPts,
        maxEqOp<point>(),
        -greatPoint
    );

    //Reset layer points into newPoints field
    forAll(tPts, pointi)
    {
        if (tPts[pointi] != -greatPoint)
        {
            newPoints[pointi] = tPts[pointi];
        }
    }

    tPts = (newPoints-mesh.points());

    syncTools::syncPointList
    (
        mesh,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = mesh.points() + tPts;
}


void Foam::layerManipulate::adjustNumLayers
(
    const indirectPrimitivePatch& pp,
    const List<Tuple2<word,scalar>>& layerMethod,
    const scalarField& initExpansion,
    const pointField& pointOrigin,
    const pointField& outerPos,
    const labelList& totalNumLayers,
    labelList& adjustedNLayers
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& meshPoints = pp.meshPoints();
    scalarField minTargetExp(meshPoints.size(), GREAT);
    scalarField maxTargetExp(meshPoints.size(), -GREAT);

    const List<Tuple2<scalar,scalar>>& targetExpansions =
       layerParams_.targetExpansions();

    forAll(pp, i)
    {
        label meshfacei = pp.addressing()[i];
        label patchi = patches.whichPatch(meshfacei);
        const face& f = pp.localFaces()[i];
        const Tuple2<scalar,scalar>& targetExp = targetExpansions[patchi];
        scalar minExp = targetExp.first();
        scalar maxExp = targetExp.second();
        forAll(f,fp)
        {
            label pointi = f[fp];
            minTargetExp[pointi] = min(minTargetExp[pointi],minExp);
            maxTargetExp[pointi] = max(maxTargetExp[pointi],maxExp);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        minTargetExp,
        minEqOp<scalar>(),
        -GREAT         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        maxTargetExp,
        maxEqOp<scalar>(),
        GREAT         // null value
    );

    forAll(meshPoints, ptI)
    {
        label meshPointI = meshPoints[ptI];
        if (layerMethod[ptI].first() == "fch")
        {
            scalar eR = initExpansion[ptI];
            scalar expMin = minTargetExp[ptI];
            scalar expMax = maxTargetExp[ptI];
            scalar expTarget = -1;
            if (eR < expMin)
            {
                expTarget = expMin;
            }
            if (eR > expMax)
            {
                expTarget = expMax;
            }
            label numLayers = totalNumLayers[meshPointI];
            scalar firstCellHeight = layerMethod[ptI].second();
            vector layerVec = outerPos[meshPointI]-pointOrigin[meshPointI];
            scalar layerHeight = mag(layerVec);
            label wantedLayers = numLayers;
            if (wantedLayers > 1 && expTarget > -1)
            {
               if (expTarget < 1.0 - SMALL || expTarget > 1.0 + SMALL)
               {
                   scalar a1 = 1.
                       - (layerHeight/firstCellHeight)*(1. - expTarget);
                   wantedLayers = max(label(2),log(a1)/log(expTarget));
               }
               else
               {
                   wantedLayers = max(label(2),layerHeight/firstCellHeight);
               }
            }
            adjustedNLayers[ptI] = wantedLayers;
        }
    }

    return;
}


void Foam::layerManipulate::reProjectOuter
(
    const labelList& ftrPointOrigin,
    const labelList& layerPointType,
    const scalarField& maxLayerThickness,
    const boolList& boundaryEdges,
    const boolList& stackEdges,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    const pointField& newPoints,
    pointField& outerPos
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    boolList pointVisited(mesh.nPoints(), false);
    forAll(mesh.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            pointVisited[pointi] = true;

            if (ftrPointOrigin[pointi] == 2)
            {
                point startPt = pointOrigin[pointi];
                scalar layerHeight = mag(newPoints[pointi]-startPt);
                point endPt = startPt
                    + 2.*layerHeight*pointSurfNormal[pointi];
                pointHit lHit = linePointRef(startPt, endPt).nearestDist
                    (newPoints[pointi]);
                if (lHit.hit())
                {
                    outerPos[pointi] =
                        0.5*(lHit.hitPoint()+newPoints[pointi]);
                }
                else
                {
                    outerPos[pointi] = newPoints[pointi];
                }
            }
            else
            {
                outerPos[pointi] = newPoints[pointi];
            }
        }
    }

    limitLayerHeight
    (
        ftrPointOrigin,
        layerPointType,
        maxLayerThickness,
        pointSurfNormal,
        pointOrigin,
        outerPos
    );

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh,
            outerPos,
            maxMagSqrEqOp<point>(),
            vector::zero
         );

        syncTools::syncPointList
        (
            mesh,
            pointVisited,
            orEqOp<bool>(),
            false
         );

        forAll(mesh.edges(), edgei)
        {
            if (!boundaryEdges[edgei] && stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                if
                (
                    pointVisited[e[0]] && !pointVisited[e[1]]
                )
                {
                    pointVisited[e[1]] = true;
                    outerPos[e[1]] = outerPos[e[0]];
                    nSet++;
                }
                else if
                (
                    pointVisited[e[1]] && !pointVisited[e[0]]
                )
                {
                    pointVisited[e[0]] = true;
                    outerPos[e[0]] = outerPos[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }
}


void Foam::layerManipulate::reSnap
(
    const indirectPrimitivePatch& grownUpSnapPP,
    const autoPtr<indexedOctree<treeDataEdge>>& grownUpFeatureMeshPtr,
    const indirectPrimitivePatch& fzonePP,
    const autoPtr<searchableSurfaces>& grownUpGeometryPtr,
    const autoPtr<searchableSurfaces>& grownUpZoneGeometryPtr,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const vectorField& newFaceAreas,
    const pointField& newFaceCentres,
    pointField& tPts,
    pointField& newPoints
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const scalar edge0Length = meshRefiner_.meshCutter().level0EdgeLength();

    //Re-snap if any patches selected
    if (returnReduce(grownUpSnapPP.size(), sumOp<label>()) != 0)
    {
        const labelList& grownUpSnapPPMeshPts = grownUpSnapPP.meshPoints();

        labelList snapSurf = identity(grownUpGeometryPtr().size());
        pointField snapPts(grownUpSnapPPMeshPts.size());
        scalarField snapDist(grownUpSnapPPMeshPts.size());
        forAll(snapPts, ptI)
        {
            label meshPointI = grownUpSnapPPMeshPts[ptI];
            snapPts[ptI] = newPoints[meshPointI];
            snapDist[ptI] =
                sqr(edge0Length / (1<<pointLevel[meshPointI]));
        }

        List<pointIndexHit> info;
        grownUpGeometryPtr().findNearest
        (
            snapPts,
            snapDist,
            snapSurf,
            info
        );

        bool snapToFtr(grownUpFeatureMeshPtr.valid() ? true : false);

        forAll(snapPts, ptI)
        {
            label meshPointI = grownUpSnapPPMeshPts[ptI];
            if (snapToFtr && boundaryPoints[meshPointI] == 2)
            {
                point sample = newPoints[meshPointI];
                scalar ftrDistSqr = sqr
                (
                    edge0Length / (1<<pointLevel[meshPointI])
                );
                pointIndexHit ftrInfo = grownUpFeatureMeshPtr().findNearest
                (
                    sample,
                    ftrDistSqr
                );
                if (ftrInfo.hit())
                {
                    tPts[meshPointI] += ftrInfo.hitPoint()
                        - newPoints[meshPointI];
                }
                else if (info[ptI].hit())
                {
                    tPts[meshPointI] += info[ptI].hitPoint()
                        - newPoints[meshPointI];
                }
            }
            else if (info[ptI].hit())
            {
                if
                (
                    boundaryPoints[meshPointI] == 1
                    || boundaryPoints[meshPointI] == 2
                )
                {
                    tPts[meshPointI] += info[ptI].hitPoint()
                        - newPoints[meshPointI];
                }
            }
        }
    }

    if
    (
        grownUpZoneGeometryPtr().size()
        && returnReduce(fzonePP.meshPoints().size(), sumOp<label>()) != 0
    )
    {
        const labelList& zoneMeshPts = fzonePP.meshPoints();
        pointField snapPts(zoneMeshPts.size());
        scalarField snapDist(zoneMeshPts.size());

        labelList snapSurf = identity(grownUpZoneGeometryPtr().size());
        pointField ppNewZonePts(zoneMeshPts.size(), vector::zero);
        scalarField ppNewZoneWeights(zoneMeshPts.size(), Zero);

        //locally smooth zone points
        forAll(snapPts, ptI)
        {
            if
            (
                stationaryPts[zoneMeshPts[ptI]] > -1
                || boundaryPoints[zoneMeshPts[ptI]] > -1
             )
            {
                continue;
            }

            const labelList& pFaces = fzonePP.pointFaces()[ptI];
            forAll(pFaces, pfI)
            {
                label facei = fzonePP.addressing()[pFaces[pfI]];
                scalar fArea = mag(newFaceAreas[facei]);
                ppNewZonePts[ptI] += fArea*newFaceCentres[facei];
                ppNewZoneWeights[ptI] += fArea;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            zoneMeshPts,
            ppNewZonePts,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            zoneMeshPts,
            ppNewZoneWeights,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
        );

        forAll(snapPts, ptI)
        {
            label meshPointI = zoneMeshPts[ptI];
            point avePt = newPoints[meshPointI];
            if (ppNewZoneWeights[ptI] > VSMALL)
            {
                avePt = ppNewZonePts[ptI]/ppNewZoneWeights[ptI];
            }

            snapPts[ptI] = avePt;
            snapDist[ptI] =
                sqr(edge0Length / (1<<pointLevel[meshPointI]));
        }

        List<pointIndexHit> info;
        grownUpZoneGeometryPtr().findNearest
        (
            snapPts,
            snapDist,
            snapSurf,
            info
         );

        forAll(snapPts, ptI)
        {
            if (info[ptI].hit())
            {
                label meshPointI = zoneMeshPts[ptI];
                if
                (
                    (
                       boundaryPoints[meshPointI] == -1
                       || boundaryPoints[meshPointI] == 1
                    )
                    &&
                    (
                        stationaryPts[meshPointI] < 0
                        || stationaryPts[meshPointI] == 3
                    )
                 )
                {
                    tPts[meshPointI] += info[ptI].hitPoint()
                        - newPoints[meshPointI];
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = newPoints + tPts;
}


void Foam::layerManipulate::correctConcaveCells
(
    const labelList& layerPointType,
    const labelList& layerFaceType,
    pointField& newFaceCentres,
    pointField& newCellCentres,
    vectorField& newFaceAreas,
    scalarField& newCellVolumes,
    pointField& tPts,
    pointField& newPoints
)
{
    tPts = vector::zero;

    fvMesh& mesh = meshRefiner_.mesh();
    const labelIOList& layerCells = meshRefiner_.layerCells()();

    updateGeometry
    (
        newPoints,
        newFaceCentres,
        newFaceAreas,
        newCellCentres,
        newCellVolumes
    );

    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);
    labelList checkFaces(identity(mesh.nFaces()));
    Foam::polyMeshGeometry::checkFacePyramids
    (
        false,
        -VSMALL,
        mesh,
        newFaceCentres,
        newCellCentres,
        newFaceAreas,
        newPoints,
        checkFaces,
        List<labelPair>(0),
        &errorFaces
    );

    if (returnReduce(errorFaces.size(), sumOp<label>()) == 0)
    {
       return;
    }

    labelList nMoved(mesh.nPoints(), 0);
    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] < 0)
        {
            const cell& c = mesh.cells()[celli];
            DynamicList<label> stackFaces(c.size());
            bool cellHasErrors = false;
            DynamicList<label> cellErrorFaces(c.size());
            forAll(c, cFI)
            {
                label facei = c[cFI];
                if (layerFaceType[facei] == 3)
                {
                    stackFaces.append(facei);
                }
                if (errorFaces.found(facei))
                {
                    cellHasErrors = true;
                }
            }

            if (stackFaces.size() > 1 && cellHasErrors)
            {
                labelHashSet stackFaceSet(stackFaces);
                const labelList& cEdges = mesh.cellEdges()[celli];
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    const edge& e = mesh.edges()[edgei];
                    if
                    (
                        layerPointType[e[0]] == 2
                        && layerPointType[e[1]] == 2
                    )
                    {
                        const labelList& eFaces =
                            mesh.edgeFaces()[edgei];
                        DynamicList<label> cellEdgeFaces(2);
                        DynamicList<label> cellErrorFaces(2);
                        forAll(eFaces, eFI)
                        {
                            label facei = eFaces[eFI];
                            if
                            (
                                stackFaceSet.found(facei)
                            )
                            {
                                cellEdgeFaces.append(facei);
                                if (errorFaces.found(facei))
                                {
                                   cellErrorFaces.append(facei);
                                }
                            }
                        }
                        if
                        (
                            cellEdgeFaces.size() == 2
                            && cellErrorFaces.size() == 1
                        )
                        {
                            pointField fCentres(2);
                            pointField fNormals(2);
                            scalarField fAreas(2);
                            forAll(cellEdgeFaces, cEFI)
                            {
                                label facei = cellEdgeFaces[cEFI];
                                const face& f = mesh.faces()[facei];
                                label own = mesh.faceOwner()[facei];
                                fNormals[cEFI] = f.areaNormal(newPoints);
                                if (own == celli)
                                {
                                    fNormals[cEFI] = -fNormals[cEFI];
                                }
                                fCentres[cEFI] = f.centre(newPoints);
                                fAreas[cEFI] = mag(fNormals[cEFI]);
                                fNormals[cEFI] /= (fAreas[cEFI]+VSMALL);
                            }
                            scalar sumArea = fAreas[0]+fAreas[1];

                            vector eNorm
                            (
                                fAreas[0]*fNormals[0]
                                +fAreas[1]*fNormals[1]
                            );
                            eNorm /= (sumArea + VSMALL);

                            if (mag(eNorm) < VSMALL)
                            {
                                continue;
                            }

                            point eC = e.centre(newPoints);
                            vector fC2eC = fCentres[0] - eC;
                            //check if concave cell
                            if ((fC2eC&eNorm) < 0)
                            {
                               label errorFace = cellErrorFaces[0];
                               label otherFace
                               (
                                  cellEdgeFaces[0] == errorFace
                                  ?  cellEdgeFaces[1] : cellEdgeFaces[0]
                               );
                               const face& f = mesh.faces()[errorFace];
                               point otherFC = newFaceCentres[otherFace];
                               point otherNorm = newFaceAreas[otherFace];
                               plane fPlane(otherFC, otherNorm);
                               forAll(f,fp)
                               {
                                  label pointi = f[fp];
                                  if (pointi != e[0] && pointi != e[1])
                                  {
                                     vector disp = fPlane.nearestPoint
                                     (
                                        newPoints[pointi]
                                     ) - newPoints[pointi];
                                     tPts[pointi] += disp;
                                     nMoved[pointi]++;
                                  }
                               }
                            }
                        }
                    }
                }
            }
        }
    }

    forAll(nMoved, pointi)
    {
        if (nMoved[pointi] > 0)
        {
            point avePt = newPoints[pointi]
                + (tPts[pointi]/nMoved[pointi]);
            newPoints[pointi] = avePt;
        }
    }
}


void Foam::layerManipulate::resetPoints
(
    const label iter,
    const meshControl& controller,
    const labelList& layerPointType,
    const labelList& layerFaceType,
    const boolList& layerEdges,
    const scalarField& medialDist,
    const scalarField& medialLayerHeight,
    const pointField& pointOrigin,
    const pointField& pointSurfNormal,
    const pointField& tPts,
    const pointField& outerPos,
    scalarField& newCellVolumes,
    vectorField& newFaceAreas,
    pointField& newFaceCentres,
    pointField& newCellCentres,
    pointField& newPoints,
    labelList& resetPts,
    bool correctErrors
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const scalar edge0Length = meshRefiner_.meshCutter().level0EdgeLength();

    if (iter > 0)
    {
        forAll(mesh.points(), pointi)
        {
            if (layerPointType[pointi] == 2)
            {
                vector toSurf = newPoints[pointi] - pointOrigin[pointi];
                scalar dist0 = mag(toSurf);

                if
                (
                    dist0 > 0.8*medialDist[pointi]
                    && dist0 < 1.1*medialLayerHeight[pointi]
                 )
                {
                    newPoints[pointi] = tPts[pointi];
                }
            }
        }
    }

    forAll(mesh.points(), pointi)
    {
        if (resetPts[pointi] == 1)
        {
            newPoints[pointi] = tPts[pointi];
            const labelList& pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (!layerEdges[edgei])
                {
                    edge e = mesh.edges()[edgei];
                    label otherPt = (e[0] ==  pointi ? e[1] : e[0]);
                    if (layerPointType[otherPt] >= 2)
                    {
                        newPoints[otherPt] =
                            0.5*(tPts[otherPt]+newPoints[otherPt]);
                    }
                }
            }
        }
    }

    //If errors introduced reset outer layer point
    if (correctErrors)
    {
        scalar maxOrtho = layerParams_.dualOrtho();

        const labelIOList& layerCells = meshRefiner_.layerCells()();

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);
        labelList checkFaces(identity(mesh.nFaces()));

        Foam::polyMeshGeometry::checkFacePyramids
        (
            false,
            -VSMALL,
            mesh,
            newFaceCentres,
            newCellCentres,
            newFaceAreas,
            newPoints,
            checkFaces,
            List<labelPair>(0),
            &errorFaces
         );

        if (maxOrtho < 180.0-SMALL)
        {
            Foam::polyMeshGeometry::checkFaceDotProduct
            (
                false,
                maxOrtho,
                scalar(180),
                mesh,
                newFaceCentres,
                newCellCentres,
                newFaceAreas,
                checkFaces,
                List<labelPair>(0),
                &errorFaces
            );
        }

        boolList errorCells(mesh.nCells(), false);
        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label facei = iter.key();
            label own = mesh.faceOwner()[facei];
            if (layerCells[own] > -1)
            {
                errorCells[own] = true;
            }

            if (mesh.isInternalFace(facei))
            {
                label nei = mesh.faceNeighbour()[facei];
                if (layerCells[nei] > -1)
                {
                    errorCells[nei] = true;
                }
            }
        }

        boolList errorPts(mesh.nPoints(), false);
        forAll(errorCells, celli)
        {
            if (errorCells[celli])
            {
                const labelList& cellPts = mesh.cellPoints()[celli];
                forAll(cellPts, cPI)
                {
                    errorPts[cellPts[cPI]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            errorPts,
            orEqOp<bool>(),
            false
         );

        forAll(errorPts, pointi)
        {
            if (errorPts[pointi])
            {
                const labelList& pCells = mesh.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    if (layerCells[celli] > -1)
                    {
                        errorCells[celli] = true;
                    }
                }
            }
        }

        labelList faceLayerCells(mesh.nFaces(), -1);
        forAll(errorCells, celli)
        {
            if (errorCells[celli])
            {
                label layerCellI = layerCells[celli];
                const cell& c = mesh.cells()[celli];
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    faceLayerCells[facei] =
                        max(faceLayerCells[facei],layerCellI);
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            faceLayerCells,
            maxEqOp<label>()
        );

        while (true)
        {
            label nReset = 0;
            forAll(errorCells, celli)
            {
                if (!errorCells[celli] && layerCells[celli] > -1)
                {
                    label layerCellI = layerCells[celli];
                    const cell& c = mesh.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        if (faceLayerCells[facei] == layerCellI)
                        {
                            errorCells[celli] = true;
                            break;
                        }
                    }
                }
            }

            if (returnReduce(nReset, sumOp<label>()) == 0)
            {
                break;
            }

            forAll(errorCells, celli)
            {
                if (errorCells[celli])
                {
                    label layerCellI = layerCells[celli];
                    const cell& c = mesh.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        faceLayerCells[facei] =
                            max(faceLayerCells[facei],layerCellI);
                    }
                }
            }

            syncTools::syncFaceList
            (
                mesh,
                faceLayerCells,
                maxEqOp<label>()
             );
        }

        forAll(mesh.cells(), celli)
        {
            if (errorCells[celli] && layerCells[celli] > -1)
            {
                const labelList& cPts = mesh.cellPoints()[celli];
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    resetPts[pointi] = label(0);
                }
            }
        }

        //Check first non layer cells becoming too deformed (acute angles)
        if (controller.algorithm() == meshControl::EXTRUDE)
        {
            forAll(mesh.cells(), celli)
            {
                if (layerCells[celli] < 0)
                {
                    const cell& c = mesh.cells()[celli];
                    label nLayerFaces = 0;
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        const face& f = mesh.faces()[facei];

                        if (layerFaceType[facei] == 3)
                        {
                            if (f.size() == 3)
                            {
                                nLayerFaces++;
                            }
                            else
                            {
                                nLayerFaces += 2;
                            }
                        }
                    }

                    if (nLayerFaces > 2)
                    {
                        const labelList& cPts = mesh.cellPoints()[celli];
                        point newCC = newCellCentres[celli];
                        scalar cellLength = edge0Length
                            / (1<<cellLevel[celli]);
                        bool warpedCell = false;

                        forAll(cPts, cPtI)
                        {
                            label pointi = cPts[cPtI];
                            if (layerPointType[pointi] >=2)
                            {
                                point pN = pointSurfNormal[pointi];
                                scalarField vProj(cPts.size(), Zero);
                                forAll(cPts, cPtJ)
                                {
                                    label pointj = cPts[cPtJ];
                                    vector n = newPoints[pointj] - newCC;
                                    vProj[cPtJ] = (pN & n);
                                }
                                // Get normal 'span' of cell
                                scalar minVal = min(vProj);
                                scalar maxVal = max(vProj);
                                scalar dW = (maxVal - minVal);
                                if (dW < 0.2*cellLength)
                                {
                                    warpedCell = true;
                                    break;
                                }
                            }
                        }
                        if (warpedCell)
                        {
                            forAll(cPts, cPtI)
                            {
                                resetPts[cPts[cPtI]] = label(0);
                            }
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            resetPts,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(newPoints, pointi)
        {
            if (resetPts[pointi] == 0 && layerPointType[pointi] >= 2)
            {
                newPoints[pointi] = outerPos[pointi];
            }
        }
    }
}


void Foam::layerManipulate::projectOuterLayerFaces
(
    const label& iter,
    const autoPtr<indirectPrimitivePatch>& outerShellPtr,
    const labelList& layerPointType,
    const labelList& layerFaceType,
    const labelList& ftrPointOrigin,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const boolList& stackEdges,
    const PackedBoolList& isMasterFace,
    const vectorField& newFaceAreas,
    const pointField& newFaceCentres,
    const pointField& newCellCentres,
    const pointField& pointOrigin,
    const vectorField& pointSurfNormal,
    const pointField& tPts,
    pointField& newPoints,
    labelList& resetPts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();

    const labelIOList& layerCells = meshRefiner_.layerCells()();

    pointField neiNewCellCentres(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiNewCellCentres[facei-mesh.nInternalFaces()] =
            newCellCentres[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiNewCellCentres);

    scalarField weightSum(mesh.nPoints(), 0);

    forAll(mesh.faces(), facei)
    {
        if
        (
            !isMasterFace[facei] || layerFaceType[facei] != 3
        )
        {
            continue;
        }

        label patchI = patches.whichPatch(facei);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[facei];

            point nCC = vector::zero;
            if (patchI == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                nCC = newCellCentres[nei];
            }
            else
            {
                nCC = neiNewCellCentres[facei-mesh.nInternalFaces()];
            }

            point midPoint = 0.5*(newCellCentres[own]+nCC);
            scalar wt = mag(newFaceAreas[facei])+ VSMALL;
            face f = mesh.faces()[facei];
            vector pNorm = vector::zero;
            forAll(f, fp)
            {
                label pointi = f[fp];
                pNorm += pointSurfNormal[pointi];
            }
            pNorm /= f.size();
            scalar magPNorm = mag(pNorm);

            if (magPNorm > VSMALL)
            {
                pNorm /= magPNorm;
                plane fPlane(midPoint, pNorm);
                Pair<vector> n1n2;
                n1n2[0] = newFaceAreas[facei];
                n1n2[0] /= (mag(n1n2[0]) + VSMALL);
                if (layerCells[owners[facei]] > -1)
                {
                   n1n2[0] = -n1n2[0];
                }
                n1n2[1] = pNorm;

                tensor T = rotationTensor(n1n2[0], -n1n2[1]);
                point currFC = newFaceCentres[facei];
                forAll(f, fp)
                {
                    // rotate points
                    label pointi = f[fp];
                    point currPt = tPts[pointi];
                    if (stationaryPts[pointi] < 0)
                    {
                        currPt -= currFC;
                        currPt = transform(T, currPt);
                        currPt += currFC;
                        vector disp = fPlane.nearestPoint(currPt)
                            - tPts[pointi];
                        scalar fWeight = (1./wt);
                        weightSum[pointi] += fWeight;
                        newPoints[pointi] += fWeight*disp;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        weightSum,
        plusEqOp<scalar>(),
        scalar(0)
     );

    syncTools::syncPointList
    (
        mesh,
        newPoints,
        plusEqOp<point>(),
        vector::zero
     );

    forAll(newPoints, pointi)
    {
        if
        (
            layerPointType[pointi] >= 2 && weightSum[pointi] > VSMALL
        )
        {
            point updatedPt = tPts[pointi]
                + 0.5*(newPoints[pointi]/weightSum[pointi]);
            newPoints[pointi] = updatedPt;
        }
        else
        {
            newPoints[pointi] = tPts[pointi];
        }
    }

    //Check for inverted edges dues to projection
    if (outerShellPtr.valid())
    {
        const indirectPrimitivePatch& outerShell = outerShellPtr();

        const labelList& outerPts =  outerShell.meshPoints();
        const edgeList& outerEdges = outerShell.edges();

        if ((iter % 2) == 0)
        {
            const labelList outerMeshEdges
            (
                outerShell.meshEdges
                (
                    mesh.edges(),
                    mesh.pointEdges()
                )
            );

            autoPtr<boolList> flipMap
            (
                new boolList(outerShell.size(), false)
            );
            forAll(outerShell, i)
            {
                label facei = outerShell.addressing()[i];
                label ownLayerCell =
                   layerCells[mesh.faceOwner()[facei]];
                if (ownLayerCell > -1)
                {
                    flipMap()[i] = true;
                }
            }

            edgeClassification eOuterClass
            (
                mesh,
                mesh.points(),
                outerShell,
                outerMeshEdges,
                0.707,
                0.707,
                -GREAT,
                flipMap
            );
            const List<Tuple2<edgeClassification::edgeType,scalar>>&
                eOuterType = eOuterClass.edgeTypes();
            const scalar edge0Length =
                meshRefiner_.meshCutter().level0EdgeLength();
            const labelList& pointLevel  =
                meshRefiner_.meshCutter().pointLevel();

            // Smooth outer shell
            weightSum = 0.0;
            pointField smoothedPts(mesh.nPoints(), vector::zero);
            pointField outerFC(outerShell.size(), vector::zero);
            scalarField outerFA(outerShell.size(), 0.0);
            boolList smoothPts(outerPts.size(), false);
            forAll(outerShell, facei)
            {
                label meshfacei = outerShell.addressing()[facei];
                const face& f = mesh.faces()[meshfacei];
                outerFA[facei] = f.mag(newPoints);
                outerFC[facei] = f.centre(newPoints);

                const face& lf = outerShell.localFaces()[facei];

                // Check for points to smooth
                forAll(lf,fp)
                {
                    label nextFp = lf.fcIndex(fp);
                    label prevFp = lf.rcIndex(fp);
                    label currPt = outerPts[lf[fp]];
                    label nextPt = outerPts[lf[nextFp]];
                    label prevPt = outerPts[lf[prevFp]];
                    vector nextEdgeVec = newPoints[nextPt] - newPoints[currPt];
                    vector prevEdgeVec = newPoints[prevPt] - newPoints[currPt];
                    scalar nextEdgeMag = mag(nextEdgeVec);
                    scalar prevEdgeMag = mag(prevEdgeVec);
                    label maxLevel = max
                    (
                        pointLevel[currPt],
                        pointLevel[nextPt]
                    );

                    label lei =  meshTools::findEdge
                    (
                        outerEdges,
                        outerShell.pointEdges()[lf[fp]],
                        lf[fp],
                        lf[nextFp]
                    );

                    if (lei != -1)
                    {
                        scalar undistortedLen = edge0Length / pow(2., maxLevel);
                        if (nextEdgeMag < 0.02*undistortedLen)
                        {
                            if (ftrPointOrigin[currPt] == -1)
                            {
                                smoothPts[lf[fp]] = true;
                            }
                            if (ftrPointOrigin[nextPt] == -1)
                            {
                                smoothPts[lf[nextFp]] = true;
                            }
                        }
                        else if
                        (
                            (
                                eOuterType[lei].first()
                                == edgeClassification::CONCAVE
                            )
                            ||
                            (
                                eOuterType[lei].first()
                                == edgeClassification::CONVEX
                            )
                        )
                        {
                            if (ftrPointOrigin[currPt] == -1)
                            {
                                smoothPts[lf[fp]] = true;
                            }
                            if (ftrPointOrigin[nextPt] == -1)
                            {
                                smoothPts[lf[nextFp]] = true;
                            }
                        }
                    }

                    if (nextEdgeMag > VSMALL && prevEdgeMag > VSMALL)
                    {
                        nextEdgeVec /= nextEdgeMag;
                        prevEdgeVec /= prevEdgeMag;
                        scalar dotProd = (nextEdgeVec & prevEdgeVec);
                        if (dotProd  > 0.939)
                        {
                            smoothPts[lf[fp]] = true;
                        }
                        else if (dotProd  < -0.939)
                        {
                            smoothPts[lf[prevFp]] = true;
                            smoothPts[lf[nextFp]] = true;
                        }
                    }
                    else
                    {
                        smoothPts[lf[fp]] = true;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                outerPts,
                smoothPts,
                orEqOp<bool>(),
                false
            );

            forAll(outerPts, pti)
            {
                label meshpointi = outerPts[pti];
                if
                (
                    smoothPts[pti]
                    && stationaryPts[meshpointi] < 0
                    && boundaryPoints[meshpointi] == -1
                )
                {
                    const labelList& pFaces = outerShell.pointFaces()[pti];
                    forAll(pFaces, pFI)
                    {
                        label facei = pFaces[pFI];
                        label meshfacei = outerShell.addressing()[facei];
                        label faceLevel =
                            meshRefiner_.meshCutter().faceLevel(meshfacei);
                        scalar len = edge0Length / pow(2., faceLevel);
                        scalar fWeight = sqrt(outerFA[facei])/len;
                        weightSum[meshpointi] += fWeight;
                        smoothedPts[meshpointi] += fWeight*outerFC[facei] ;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                weightSum,
                plusEqOp<scalar>(),
                scalar(0)
            );

            syncTools::syncPointList
            (
                mesh,
                smoothedPts,
                plusEqOp<vector>(),
                vector::zero
            );

            forAll(mesh.points(), pointi)
            {
                if (weightSum[pointi] > VSMALL)
                {
                    point avePt = smoothedPts[pointi] / weightSum[pointi];
                    newPoints[pointi] = 0.5*(newPoints[pointi]+avePt);
                }
            }
        }

        forAll(outerPts, pti)
        {
            label pointi = outerPts[pti];

            const labelList& pEdges = outerShell.pointEdges()[pti];
            forAll(pEdges, pEI)
            {
                const edge& e = outerEdges[pEdges[pEI]];
                label otherpointi = outerPts[e.otherVertex(pti)];

                if
                (
                    ftrPointOrigin[pointi] != 0
                    && ftrPointOrigin[otherpointi] != 0
                 )
                {
                    continue;
                }

                vector eVec = newPoints[otherpointi]
                    - newPoints[pointi];
                scalar eVecMag  = mag(eVec);

                vector eVecOrig = tPts[otherpointi]
                    - tPts[pointi];
                scalar eVecOrigMag  = mag(eVecOrig);

                vector fVec = pointOrigin[otherpointi]
                    - pointOrigin[pointi];
                scalar fVecMag  = mag(fVec);

                if (eVecMag < VSMALL || fVecMag < VSMALL)
                {
                    resetPts[pointi] = label(0);
                    resetPts[otherpointi] = label(0);
                }
                else
                {
                    eVec /= eVecMag;
                    fVec /= fVecMag;

                    if
                    (
                        (
                            ( eVecMag < eVecOrigMag)
                            && (eVecMag < 0.25*fVecMag)
                        )
                        || (fVec & eVec) < 0
                    )
                    {
                        resetPts[pointi] = label(0);
                        resetPts[otherpointi] = label(0);
                    }
                }
            }
        }
    }

    forAll(newPoints, pointi)
    {
       if (layerPointType[pointi] >= 2)
       {
            labelList pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];

                if (stackEdges[edgei])
                {
                    edge e = mesh.edges()[edgei];
                    label otherPt = (e[0] ==  pointi ? e[1] : e[0]);

                    vector toSurf = newPoints[pointi] - pointOrigin[otherPt];
                    scalar toSurfMag = mag(toSurf);
                    toSurf /= (toSurfMag + VSMALL);

                    vector sNorm = pointSurfNormal[otherPt];
                    sNorm /= (mag(sNorm) + VSMALL);

                    vector toSurfOrig = tPts[pointi]
                        - pointOrigin[otherPt];
                    scalar toSurfOrigMag = mag(toSurfOrig);
                    toSurfOrig /= (toSurfOrigMag  + VSMALL);

                    scalar dP1 = (toSurf&sNorm);
                    scalar dP2 = (toSurfOrig&sNorm);

                    if
                    (
                       dP1 < dP2 && dP1 < 0.866
                       && ftrPointOrigin[otherPt] != -1
                    )
                    {
                        resetPts[pointi] = label(1);
                        break;
                    }
                }
            }
       }
    }

    syncTools::syncPointList
    (
        mesh,
        resetPts,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(newPoints, pointi)
    {
        if (resetPts[pointi] == 0)
        {
            newPoints[pointi] = tPts[pointi];
        }
    }
}


void Foam::layerManipulate::limitLayerHeight
(
    const labelList& ftrPointOrigin,
    const labelList& layerPointType,
    const scalarField& maxLayerThickness,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    pointField& pts
)
{
    forAll(layerPointType, pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            vector layerVec = pts[pointi]
                - pointOrigin[pointi];
            scalar lheight = mag(layerVec);
            scalar mlh = maxLayerThickness[pointi];

            if (lheight > mlh && mlh > 0)
            {
               if (ftrPointOrigin[pointi] == -1)
                {
                    vector sNorm = pointSurfNormal[pointi];
                    sNorm /= (mag(sNorm) + VSMALL);
                    plane fPlane(pointOrigin[pointi], sNorm);
                    point hPt =
                        fPlane.nearestPoint(pts[pointi]);
                    scalar nheight = mag(pts[pointi]-hPt);
                    if (nheight > mlh)
                    {
                        scalar correction(mlh/nheight);
                        pts[pointi] = pointOrigin[pointi]
                            + correction*layerVec;
                    }
                }
                else
                {
                    scalar correction(mlh/lheight);
                    pts[pointi]  = pointOrigin[pointi]
                        + correction*layerVec;
                }
            }
        }
    }
}


void Foam::layerManipulate::fitLayerPointStack()
{
    Info<<"Fitting Layer Stacks"<<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    //Create a unique tag for the points based on the labels of layerCells
    labelListList pointTagList = createLayerPointTagList();
    pointField state = mesh.points();

    boolList needForSynchPoints(mesh.nPoints(), false);
    boolList alreadyAddedToStack(mesh.nPoints(), false);
    List<layerPointStack> layerPointStackList(mesh.nPoints());


    labelList totalNuOfLayers(mesh.nPoints(), -1);
    forAll(pp.meshPoints(), pI)
    {
        const labelList& pointFaces = pp.pointFaces()[pI];
        forAll(pointFaces, fI)
        {
            const label& fL = pp.addressing()[fI];
            const label& patch = mesh.boundaryMesh().whichPatch(fL);
            totalNuOfLayers[pI] =
                max(totalNuOfLayers[pI], (layerParams_.numLayers()[patch] + 1));
        }
    }

    forAll(pp.meshPoints(), pI)
    {
        const label& pL = pp.meshPoints()[pI];
        layerPointStack lPS = buildLayerStack
        (
            pI,
            pointTagList,
            pL,
            alreadyAddedToStack,
            totalNuOfLayers
        );
        const label& size = lPS.pointStackLabels().size();
        if (size<lPS.layerNu())
        {
            const label& lastLabel = lPS.pointStackLabels()[size-1];
            needForSynchPoints[lastLabel] = true;
        }
        forAll(lPS.pointStackLabels(), i)
        {
            const label& pL = lPS.pointStackLabels()[i];
            layerPointStackList[pL] = lPS;
            layerPointStackList[pL].setAddress(pL);
        }
        const label& baseLabel = lPS.pointStackLabels()[0];
        for (int i=1; i<size;i++)
        {
            totalNuOfLayers[i] = totalNuOfLayers[baseLabel];
        }
    }


    bool fullySynced = false;
    label numberOfSyncs = 0;
    while (!fullySynced)
    {
        label nuOfSyncPointsInit = 0;
        forAll(mesh.points(),pI)
        {
            if (needForSynchPoints[pI])
            {
                nuOfSyncPointsInit++;
            }
        }
        reduce(nuOfSyncPointsInit, sumOp<label>());

        syncTools::syncPointList
        (
            mesh,
            needForSynchPoints,
            maxEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh,
            totalNuOfLayers,
            maxEqOp<label>(),
            label(0)
        );

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& polyP = mesh.boundaryMesh()[patchI];
            if (polyP.coupled())
            {
                forAll(polyP.meshPoints(), pI)
                {
                    label pL = polyP.meshPoints()[pI];
                    if (needForSynchPoints[pL] && !alreadyAddedToStack[pL])
                    {
                        layerPointStack lPS = buildLayerStack
                        (
                            pI,
                            pointTagList,
                            pL,
                            alreadyAddedToStack,
                            totalNuOfLayers
                        );

                        forAll(lPS.pointStackLabels(), i)
                        {
                            const label& pLabel = lPS.pointStackLabels()[i];
                            layerPointStackList[pLabel] = lPS;
                            layerPointStackList[pLabel].setAddress(pL);
                        }
                        const label& size = lPS.pointStackLabels().size();
                        if (size<lPS.layerNu())
                        {
                            const label&
                                lastLabel =lPS.pointStackLabels()[size-1];
                            needForSynchPoints[lastLabel] = true;
                        }
                        const label& baseLabel = lPS.pointStackLabels()[0];
                        for (int i=1; i<size;i++)
                        {
                            totalNuOfLayers[i] = totalNuOfLayers[baseLabel];
                        }
                    }
                }
            }
        }

        label nuOfSyncPointsFinal = 0;
        forAll(mesh.points(),pI)
        {
            if (needForSynchPoints[pI])
            {
                nuOfSyncPointsFinal++;
            }
        }
        reduce(nuOfSyncPointsFinal, sumOp<label>());

        label newSyncPoints = nuOfSyncPointsFinal - nuOfSyncPointsInit;
        reduce(newSyncPoints, sumOp<label>());

        if (newSyncPoints==0)
        {
            fullySynced = true;
        }
        numberOfSyncs++;
    }

    syncTools::syncPointList
    (
        mesh,
        alreadyAddedToStack,
        maxEqOp<bool>(),
        false
    );

    List<List<point>> stackListPoints(mesh.nPoints());
    for (int i=0;i<numberOfSyncs;i++)
    {
        forAll(mesh.points(), pI)
        {
            const layerPointStack& lPS = layerPointStackList[pI];
            stackListPoints[pI] = lPS.stackPoints();
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stackListPoints,
        mergePointStacks(),
        List<point>(0)
    );

    forAll(stackListPoints, pI)
    {
        layerPointStack& lPS = layerPointStackList[pI];
        lPS.substitute(stackListPoints[pI]);
    }

    forAll(layerPointStackList, stackI)
    {
        layerPointStack& lPS = layerPointStackList[stackI];
        lPS.merge(layerPointStackList);
    }

    labelList stackIndexing(mesh.nPoints(), -1);

    forAll(layerPointStackList, stackI)
    {
        const layerPointStack& lPS = layerPointStackList[stackI];
        stackIndexing[stackI] = lPS.getIndex(state[stackI]);
    }

    forAll(layerPointStackList, stackI)
    {
        const layerPointStack& lPS = layerPointStackList[stackI];
        if (lPS.size()>0 && stackIndexing[stackI]==-1)
        {
            Pout<<" COULD NOT FIND POINT, CHECK TOLERANCE"<<endl;
            Pout<<" SEARCHING FOR "<<state[stackI]<<endl;
            Pout<<lPS.stackPoints()<<endl;
        }
    }
    forAll(layerPointStackList, stackI)
    {
        if (stackIndexing[stackI]>-1)
        {
            const layerPointStack& lPS = layerPointStackList[stackI];
            label size = lPS.size();
            pointField pF(size, vector::zero);
            forAll(pF, i)
            {
                pF[i] = lPS.stackPoints()[i];
            }
            leastSquaresCurveFit fit(pF);
            pointField fitPoints = fit.getNewPoints();
            state[stackI] = fitPoints[stackIndexing[stackI]];
        }
    }

    mesh.movePoints(state);
}

Foam::labelListList Foam::layerManipulate::createLayerPointTagList()
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelIOList& layerCells = meshRefiner_.layerCells()();

    labelListList  pointTagList(mesh.points().size());

    forAll(mesh.points(), pI)
    {
        const labelList& pointCells = mesh.pointCells()[pI];
        DynamicList<label> pointTag;
        forAll(pointCells, cI)
        {
            const label& cL = pointCells[cI];
            if (layerCells[cL] > -1)
            {
                if (pointTag.size() == 0)
                {
                    pointTag.append(layerCells[cL]);
                }
                else
                {
                    label counter = 0;
                    forAll(pointTag, elem)
                    {
                        if (pointTag[elem] != layerCells[cL])
                        {
                            counter++;
                        }
                    }
                    if (counter == pointTag.size())
                    {
                        pointTag.append(layerCells[cL]);
                    }
                }
            }
        }
        SortableList<label> sorted(pointTag);
        pointTagList[pI] = sorted;
    }

    labelList zeroSizeList(0);

    syncTools::syncPointList
    (
        mesh,
        pointTagList,
        appendUniqueElements(),
        zeroSizeList
    );

    forAll(mesh.points(), pI)
    {
        SortableList<label> sorted(pointTagList[pI]);
        pointTagList[pI] = sorted;
    }

    return pointTagList;
}

Foam::layerPointStack Foam::layerManipulate::buildLayerStack
(
    const label& basePoint,
    const labelListList& pointTagList,
    const label& ppL,
    boolList& addedToStack,
    labelList& totalNuOfLayers
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    label pL = ppL;
    label prevL = -1;
//    label counter = 1;

    DynamicList<vector> pointStack(1, mesh.points()[pL]);
    DynamicList<label> pointStackLabels(1, pL);

    addedToStack[basePoint] = true;

    bool finished = false;

    while (!finished)
    {
        finished = true;
        const labelList& pointPoints = mesh.pointPoints()[pL];
        forAll(pointPoints, pPI)
        {
            const label& pPL = pointPoints[pPI];
            const labelList& pointTagInit = pointTagList[pL];
            const labelList& pointTagNext = pointTagList[pPL];
            bool matchedTags = compareTags(pointTagInit, pointTagNext);
            if (matchedTags && pPL!=prevL)
            {
//                counter++;
                pointStack.append(mesh.points()[pPL]);
                pointStackLabels.append(pPL);
                addedToStack[pPL] = true;
                prevL = pL;
                pL = pPL;
                finished = false;
                break;
            }
        }
    }
    const label& layerNu = totalNuOfLayers[basePoint];
    layerPointStack lPS(true, pointStack, pointStackLabels, layerNu);
    return lPS;
}

bool Foam::layerManipulate::compareTags
(
    const labelList& pt1,
    const labelList& pt2
)
{
    bool matchedTags = false;
    if (pt1.size() == pt2.size())
    {
        forAll(pt1, tag)
        {
            matchedTags = true;
            if
            (
                pt2.size() > 0
             && pt1[tag]!=pt2[tag]
            )
            {
                matchedTags = false;
                break;
            }
        }
    }
    return matchedTags;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerManipulate::layerManipulate
(
    meshRefinement& meshRefiner,
    const layerParameters& layerParams
)
:
    meshRefiner_(meshRefiner),
    layerParams_(layerParams)
{}


// ************************************************************************* //
