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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "base/fvMeshGIBChangersBase.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "GIB/twoDPointGIBCorrector/twoDPointGIBCorrector.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangersBase::makeConcavePoints() const
{
    if (concavePointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList concavePoints(mesh().points().size(), false);
    boolList concaveEdges(mesh().edges().size(), false);

    const labelListList& fEdges = mesh().faceEdges();
    const labelListList& eFaces = mesh().edgeFaces();
    const vectorField& baseCf = *baseCf_;
    const pointField& baseP = *basePoints_;
    const vectorField& baseSf = *baseSf_;

    boolList boundaryEdges(mesh().edges().size(), false);
    for (label fI = mesh().nInternalFaces(); fI < mesh().nFaces(); fI++)
    {
        const labelList& fEdgesI = fEdges[fI];
        forAll(fEdgesI, eI)
        {
            label gEI = fEdgesI[eI];
            boundaryEdges[gEI] = true;
        }
    }

    List<DynamicList<label>> dbEf(mesh().edges().size());

    forAll(eFaces, eI)
    {
        if (boundaryEdges[eI])
        {
            const labelList& eFacesI = eFaces[eI];
            forAll(eFacesI, fI)
            {
                if (eFacesI[fI] >= mesh().nInternalFaces())
                {
                    dbEf[eI].append(eFacesI[fI]);
                }
            }

        }
        dbEf[eI].shrink();
    }

    forAll(boundaryEdges, eI)
    {
        if (boundaryEdges[eI])
        {
            const edge& edgeI = mesh().edges()[eI];
            const label& face1 = dbEf[eI][0];
            const label& face2 = dbEf[eI][1];
            const point& pA = baseP[edgeI.start()];
            const point& pB = baseP[edgeI.end()];

            point bC1 = baseCf[face1];
            point bC2 = baseCf[face2];

            vector triN1 = triPointRef(pA, pB, bC1).areaNormal();
            vector triN2 = triPointRef(pA, bC2, pB).areaNormal();

            if ((triN1&baseSf[face1]) < 0)
            {
                triN1 = baseSf[face1];
                triN2 = baseSf[face2];
            }
            vector AB = pB - pA;
            vector ABn = AB/mag(AB);

            vector AC1 = bC1 - pA;
            vector AC_AB1 = (AC1&ABn)*ABn;
            vector CC_AB1 = AC1 - AC_AB1;


            vector AC2 = bC2 - pA;
            vector AC_AB2 = (AC2&ABn)*ABn;
            vector CC_AB2 = AC2 - AC_AB2;

            vector CC_AB = (CC_AB1 + CC_AB2)/2;
            vector normals = (triN1 + triN2)/2;

            if ((normals&CC_AB) >0)
            {
                if (debug)
                {
                    Info<< pA <<tab << pB <<endl;
                    Info<< triN1 << tab << bC1 <<endl;
                    Info<< triN2 << tab << bC2 <<endl;
                    Info<< CC_AB << endl;
                    Info<< normals << endl;
                }

                concaveEdges[eI] = true;
                concavePoints[edgeI.start()] = true;
                concavePoints[edgeI.end()] = true;
            }
        }
    }

    concavePointsPtr_ = new boolList(concavePoints);
    concaveEdgesPtr_ = new boolList(concaveEdges);
}


void Foam::fvMeshGIBChangersBase::makeCloseBoundaryPoints() const
{
    if (closeBoundaryPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    closeBoundaryPointsPtr_ = new labelListList(mesh().points().size());
    labelListList& closeBoundaryPoints = *closeBoundaryPointsPtr_;

    List<DynamicList<label>> cBp(mesh().points().size());

    const labelListList& pointEdges = mesh().pointEdges();
    const edgeList& edges = mesh().edges();
    const labelListList& pFaces = mesh().pointFaces();

    forAll(cBp, pI)
    {
        const labelList& pEdgeI = pointEdges[pI];
        forAll(pEdgeI, eI)
        {
            label geI = pEdgeI[eI];
            const edge& edgeI = edges[geI];

            label point2 = -1;
            if (edgeI.start() == pI)
            {
                point2 = edgeI.end();
            }
            else
            {
                point2 = edgeI.start();
            }

            const labelList& pFacesI = pFaces[point2];
            forAll(pFacesI, fI)
            {
                const label& faceI = pFacesI[fI];
                if (!(faceI < mesh().nInternalFaces()))
                {
                    label patchI = mesh().boundaryMesh().whichPatch(faceI);
                    const polyPatch& poly = mesh().boundary()[patchI].patch();
                    if (isA<directPolyPatch>(poly))
                    {
                        if
                        (
                            !(
                                mesh().boundary()[patchI].coupled()
                             || isA<wedgeFvPatch>(mesh().boundary()[patchI])
                             || isA<emptyFvPatch>(mesh().boundary()[patchI])
                            )
                        )
                        {
                            bool exists(false);
                            forAll(cBp[pI], cBpI)
                            {
                                if (faceI == cBp[pI][cBpI])
                                {
                                    exists = true;
                                }
                            }
                            if (!exists)
                            {
                                cBp[pI].append(faceI);
                            }
                        }
                    }
                }
            }
        }
        cBp[pI].shrink();
    }

    forAll(closeBoundaryPoints, pI)
    {
        closeBoundaryPoints[pI] = cBp[pI];
    }
}


void Foam::fvMeshGIBChangersBase::makeMarkedBoundaryPoints() const
{
    if (markedBoundaryPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    markedBoundaryPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& markedBoundaryPoints = *markedBoundaryPointsPtr_;

    forAll(mesh().boundary(), pI)
    {
        const polyPatch& pp = mesh().boundary()[pI].patch();
        if
        (
            !(isA<indirectPolyPatch>(pp) || pp.coupled())
        )
        {
            const labelList& bPoints =
                mesh().boundary()[pI].patch().meshPoints();

            forAll(bPoints, pI)
            {
                const label& pointI = bPoints[pI];
                const labelList& pCellsI = mesh().pointCells()[pointI];
                forAll(pCellsI, pI)
                {
                    const label& cellI = pCellsI[pI];
                    const labelList& cPoints = mesh().cellPoints()[cellI];
                    forAll(cPoints, pII)
                    {
                        const label& pointI = cPoints[pII];
                        markedBoundaryPoints[pointI] = true;
                    }
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::makeNormal2DEdges() const
{
    if (normal2DEdgesPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList normal2DEdges(mesh().edges().size(), false);

    if (mesh().nGeometricD() == 2)
    {
        twoDPointGIBCorrector twoDCorr(mesh());
        const labelList& nEIndices = twoDCorr.normalEdgeIndices();
        forAll(nEIndices, eI)
        {
            const label& nEIndicesI = nEIndices[eI];
            normal2DEdges[nEIndicesI] = true;
        }
    }

    normal2DEdgesPtr_ = new boolList(normal2DEdges);
}


void Foam::fvMeshGIBChangersBase::makeBoundaryPointsRev() const
{
    if (boundaryPointsRevPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    // This is a check if the points that are on the GIB and at the boundary
    // need correction. At this point the GIB have not been updated, so the
    // check has to be made based on the faceZone.
    const boolList& interPoints = interP();
    const pointField& basePoints = *basePoints_;
    const pointField& currectPoints = mesh().points();

    boolList boundaryPointsRev = interPoints;

    forAll(boundaryPointsRev, pI)
    {
        if (boundaryPointsRev[pI])
        {
            if (basePoints[pI] == currectPoints[pI])
            {
                boundaryPointsRev[pI] = false;
            }
        }
    }

    // boolist which returns true if the edge element is internal.
    // If needed again, it should go to a new function.
    boolList internalEdge = boolList(mesh().nEdges(), true);
    const labelListList& fEdges = mesh().faceEdges();

    forAll(mesh().boundary(), pI)
    {
        const polyPatch& pp = mesh().boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if (!pp.coupled())
            {
                forAll(pp, pfI)
                {
                    label gfI = pp.start() + pfI;
                    const labelList& fEdgesI = fEdges[gfI];
                    forAll(fEdgesI, eI)
                    {
                        label geI = fEdgesI[eI];
                        internalEdge[geI] = false;
                    }
                }
            }
        }
    }

    const labelListList& pointEdges = mesh().pointEdges();
    const boolList& normal2DEdg = normal2DEdges();

    forAll(boundaryPointsRev, pI)
    {
        if (boundaryPointsRev[pI])
        {
            bool internalEdgeFound = false;
            if (!internalEdgeFound)
            {
                const labelList& pointEdgesI = pointEdges[pI];
                forAll(pointEdgesI, eI)
                {
                    label geI = pointEdgesI[eI];
                    if (internalEdge[geI] && !normal2DEdg[geI])
                    {
                        internalEdgeFound = true;
                    }
                }
            }
            if (internalEdgeFound)
            {
                boundaryPointsRev[pI] = false;
            }
        }
    }

    boundaryPointsRevPtr_ = new boolList(boundaryPointsRev);
}


void Foam::fvMeshGIBChangersBase::correctBoundaryPointsOnBaseMesh
(
    pointField& snapP
)
{
    const boolList& boundaryPointsR = boundaryPointsRev();
    const pointField& baseP = *basePoints_;
    const vectorField& baseCf = *baseCf_;

    forAll(boundaryPointsR, pI)
    {
        if (boundaryPointsR[pI])
        {
            const labelList& pFacesI = mesh().pointFaces()[pI];
            forAll(pFacesI, fI)
            {
                const label& pFI = pFacesI[fI];
                if (!(pFI < mesh().nInternalFaces()))
                {
                    label patchI = mesh().boundaryMesh().whichPatch(pFI);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchI])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchI])
                     && !isA<symmetryFvPatch>(mesh().boundary()[patchI])
                    )
                    {
                        vector dis = snapP[pI] - baseP[pI];
                        scalar dot = dis&baseCf[pFI];
                        if (dot == 0)
                        {
                            snapP[pI] = baseP[pI];
                        }
                    }
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::checkConcaveBoundaryPoints
(
    const pointField& startP,
    pointField& endP
) const
{
    const boolList& conEdges = concaveEdges();
    const labelListList& eFaces = mesh().edgeFaces();
    const vectorField& baseSf = *baseSf_;
    const edgeList& edges = mesh().edges();

    forAll(conEdges, eI)
    {
        if (conEdges[eI])
        {
            const edge& edgeI = edges[eI];
            const point& spA = startP[edgeI.start()];
            point& epA = endP[edgeI.start()];

            const point& spB = startP[edgeI.end()];
            point& epB = endP[edgeI.end()];

            vector disA = epA - spA;
            vector disB = epB - spB;

            const labelList& eFaceI = eFaces[eI];
            forAll(eFaceI, fI)
            {
                const label& faceI = eFaceI[fI];
                if (faceI >= mesh().nInternalFaces())
                {
                    label patchI = mesh().boundaryMesh().whichPatch(faceI);
                    if (!mesh().boundary()[patchI].coupled())
                    {
                        const vector& bsf = baseSf[faceI];
                        //vector bnf = (baseSf[faceI])/mag(baseSf[faceI]);

                        scalar disAn = bsf&disA;
                        scalar disBn = bsf&disB;

                        if (disAn>0)
                        {
                            //- tolerance issues
                            //epA -= (disAn*bnf);
                            epA = spA;
                        }
                        if (disBn>0)
                        {
                            //- tolerance issues
                            //epB -= (disBn*bnf);
                            epB = spB;
                        }
                    }
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::correctConstraintPatches
(
    pointField& points
) const
{
    applyTwoDCorrection(points);
    correctSymmetryPatches(points);
}


void Foam::fvMeshGIBChangersBase::correctSymmetryPatches
(
    pointField& newPoints
) const
{
    // Mark all the faces we dont want to be included
    const vectorField& baseSf = *baseSf_;
    forAll(mesh().boundary(), pI)
    {
        const polyPatch& pp = mesh().boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if
            (
                isA<symmetryFvPatch>(mesh().boundary()[pI])
            )
            {
                if (pp.size())
                {
                    const vector n = baseSf[pp.start()];

                    const labelList& pPoints = pp.meshPoints();
                    forAll(pPoints, pI)
                    {
                        label gpI = pPoints[pI];
                        newPoints[gpI] = transform(I - n*n, newPoints[gpI]);
                    }
                }
            }
        }
    }

    syncPoints(newPoints);
}


void Foam::fvMeshGIBChangersBase::applyTwoDCorrection(pointField& snapP) const
{
    if (mesh().nGeometricD() == 2)
    {
        twoDPointGIBCorrector twoDCorr(mesh());
        const pointField& basePoints = *basePoints_;
        twoDCorr.correctPointsPlus(snapP, basePoints);
        syncPoints(snapP);
    }
}


void Foam::fvMeshGIBChangersBase::applyTwoDPlanesCorrection
(
    pointField& pf1,
    pointField& pf2
)
{
    if (mesh().nGeometricD() == 2)
    {
        twoDPointGIBCorrector twoDCorr(mesh());
        twoDCorr.correctPlanes(pf1);
        twoDCorr.correctPlanes(pf2);
    }
}


// ************************************************************************* //
