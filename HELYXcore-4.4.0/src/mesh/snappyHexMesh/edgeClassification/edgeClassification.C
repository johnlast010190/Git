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
    (c) 2017-1017 Engys Ltd.
\*---------------------------------------------------------------------------*/

#include "edgeClassification/edgeClassification.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "meshTools/meshTools.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(edgeClassification, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::vectorField Foam::edgeClassification::calculatePointNormals
(
    const boolList& excludedFaces,
    const label nSmoothIter,
    const bool correctBoundaryNormals
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& meshPoints = pp_.meshPoints();
    vectorField faceNormals(pp_.size(),vector::zero);
    forAll(faceNormals, facei)
    {
        faceNormals[facei] =
            mesh_.faces()[pp_.addressing()[facei]].unitNormal(pts_);
    }

    vectorField pointNormals(pp_.nPoints(), vector::zero);
    scalarField nPointFaces(pp_.nPoints(), 0.);

    boolList externalPoints(mesh_.nPoints(), false);
    boolList externalEdges(mesh_.nEdges(), false);

    labelList nConvexEdges(pp_.nPoints(),0);
    labelList nConcaveEdges(pp_.nPoints(),0);

    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh_,
            ppMeshEdges_
         )
    );

    forAll(pp_.meshPoints(), ptI)
    {
        const labelList& pEdges = pp_.pointEdges()[ptI];
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (!isPatchMasterEdge[edgei])
            {
                continue;
            }

            if
            (
                eType_[edgei].first() == edgeClassification::CONVEX
                || eType_[edgei].first() == edgeClassification::BAFFLE
            )
            {
                nConvexEdges[ptI]++;
            }
            else if (eType_[edgei].first() == edgeClassification::CONCAVE)
            {
                nConcaveEdges[ptI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nConvexEdges,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nConcaveEdges,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(meshPoints, pointi)
    {
        label meshPointI = meshPoints[pointi];
        if (nConcaveEdges[pointi] > 3 && nConvexEdges[pointi] == 0)
        {
            const labelList& pFaces = pp_.pointFaces()[pointi];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                label meshFaceI = pp_.addressing()[facei];
                vector fN = mesh_.faces()[meshFaceI].unitNormal(pts_);
                externalPoints[meshPointI] = true;
                pointNormals[pointi] += fN;
                nPointFaces[pointi] += scalar(1.);
            }
            const labelList& pEdges = pp_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if
                (
                    eType_[edgei].first() == edgeClassification::CONCAVE
                )
                {
                    label meshedgei = ppMeshEdges_[edgei];
                    externalEdges[meshedgei] = true;
                }
            }
        }
        else if ((nConcaveEdges[pointi]+nConvexEdges[pointi]) > 0)
        {
            const labelList& pEdges = pp_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                bool concaveEdge = false;
                if (eType_[edgei].first() == edgeClassification::CONCAVE)
                {
                    concaveEdge = true;
                }

                bool convexEdge = false;
                if
                (
                    eType_[edgei].first() == edgeClassification::CONVEX
                    || eType_[edgei].first() == edgeClassification::BAFFLE
                )
                {
                    convexEdge = true;
                }

                if (concaveEdge || convexEdge)
                {
                    label meshedgei = ppMeshEdges_[edgei];
                    edge e = mesh_.edges()[meshedgei];
                    const labelList& edgeFaces = pp_.edgeFaces()[edgei];
                    forAll(edgeFaces, ppEFI)
                    {
                        label meshFaceI = pp_.addressing()[edgeFaces[ppEFI]];
                        face f = mesh_.faces()[meshFaceI];
                        vector eVec = e.vec(pts_);
                        label fp = findIndex(f, e[0]);

                        vector eNorm =
                            (eVec ^ mesh_.faces()[meshFaceI].unitNormal(pts_));

                        if (concaveEdge && f[f.fcIndex(fp)] != e[1])
                        {
                            eNorm = -eNorm;
                        }
                        else if (convexEdge && f[f.fcIndex(fp)] == e[1])
                        {
                            eNorm = -eNorm;
                        }

                        eNorm /= (mag(eNorm) + VSMALL);

                        pointNormals[pointi] += eNorm;
                        nPointFaces[pointi] += scalar(1.);
                    }

                    externalEdges[meshedgei] = true;
                    externalPoints[meshPointI] = true;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        externalEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh_,
        externalPoints,
        orEqOp<bool>(),
        false
    );

    forAll(faceNormals, facei)
    {
        if (excludedFaces[facei])
        {
            continue;
        }
        const face& f = pp_.localFaces()[facei];

        forAll(f, fp)
        {
            label meshPointI = meshPoints[f[fp]];

            if (!externalPoints[meshPointI])
            {
                scalar fA = mesh_.faces()[pp_.addressing()[facei]].mag(pts_);
                pointNormals[f[fp]] += faceNormals[facei]*fA;
                nPointFaces[f[fp]] += fA;
            }
        }
    }

    forAll(eType_, edgei)
    {
        if
        (
            eType_[edgei].first() == edgeClassification::BOUNDARY
            || eType_[edgei].first() == edgeClassification::NONMANIFOLD
        )
        {
             //Set stationary external points
            label meshedgei = ppMeshEdges_[edgei];
            edge e = mesh_.edges()[meshedgei];
            externalEdges[meshedgei] = true;
            externalPoints[e[0]] = true;
            externalPoints[e[1]] = true;
        }
    }


    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        pointNormals,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nPointFaces,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    forAll(pointNormals, i)
    {
        if (nPointFaces[i] > VSMALL)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }
    pointNormals /= (mag(pointNormals) + VSMALL);

    if (correctBoundaryNormals)
    {
        boolList markedFaces(mesh_.nFaces(), false);
        forAll(pp_, i)
        {
            markedFaces[pp_.addressing()[i]] = true;
        }

        boolList boundEdges(mesh_.nEdges(), false);
        forAll(eType_, edgei)
        {
            if (eType_[edgei].first() == edgeClassification::BOUNDARY)
            {
                label meshedgei = ppMeshEdges_[edgei];
                boundEdges[meshedgei] = true;
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            boundEdges,
            orEqOp<bool>(),
            false
        );

        labelList nNbrEdges(mesh_.nPoints(), 0);
        vectorField aveNbrDir(mesh_.nPoints(), vector::zero);
        forAll(mesh_.edges(), edgei)
        {
            if (boundEdges[edgei])
            {
                const labelList& eFaces = mesh_.edgeFaces()[edgei];
                label nbrFace = -1;
                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];
                    label patchi = patches.whichPatch(facei);
                    if
                    (
                        patchi != -1 && !patches[patchi].coupled()
                        && !markedFaces[facei]
                    )
                    {
                        nbrFace = facei;
                        break;
                    }
                }
                if (nbrFace != -1)
                {
                    const edge& e = mesh_.edges()[edgei];
                    vector eVec = e.unitVec(pts_);
                    label pt0 = e[0];
                    label pt1 = e[1];
                    const face& f = mesh_.faces()[nbrFace];
                    label start = findIndex(f, pt0);
                    if (f[f.fcIndex(start)] != pt1)
                    {
                        eVec = -eVec;
                    }
                    const vector fArea = f.unitNormal(pts_);
                    vector gUp = (eVec ^ fArea);
                    gUp /= (mag(gUp) + VSMALL);
                    aveNbrDir[pt0] += gUp;
                    aveNbrDir[pt1] += gUp;
                    nNbrEdges[pt0]++;
                    nNbrEdges[pt1]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            aveNbrDir,
            plusEqOp<vector>(),
            vector::zero        // null value
        );
        syncTools::syncPointList
        (
            mesh_,
            nNbrEdges,
            plusEqOp<label>(),
            label(0)        // null value
        );

        forAll(pointNormals, pti)
        {
            label meshpointi = meshPoints[pti];
            if (nNbrEdges[meshpointi] > 0)
            {
                point aveDir = aveNbrDir[meshpointi]
                    / nNbrEdges[meshpointi];
                aveDir /= (mag(aveDir) + VSMALL);
                pointNormals[pti] = aveDir;
            }
        }
    }

    if (debug)
    {
        pointField sPts(meshPoints.size());
        forAll(meshPoints, pti)
        {
            sPts[pti] = pts_[meshPoints[pti]];
        }
        simpleVTKWriter layerVTK
        (
            pp_.localFaces(),
            sPts
        );
        layerVTK.addPointData("pNormals", pointNormals);
        layerVTK.write("pNormUnSmooth.vtk");
    }

    if (nSmoothIter > 0)
    {
        fieldSmoother fSmoother(mesh_);

        // Precalulate (patch) master point/edge
        const PackedBoolList isPatchMasterPoint
        (
            meshRefinement::getMasterPoints
            (
                mesh_,
                meshPoints
             )
        );

        boolList cornerPoints(mesh_.nPoints(), false);
        // Smooth patch normal vectors
        fSmoother.smoothPatchNormals
        (
            nSmoothIter,
            isPatchMasterPoint,
            isPatchMasterEdge,
            externalEdges,
            externalPoints,
            cornerPoints,
            ppMeshEdges_,
            pp_,
            pointNormals
        );
    }

    if (debug)
    {
        pointField sPts(meshPoints.size());
        forAll(meshPoints, pti)
        {
            sPts[pti] = pts_[meshPoints[pti]];
        }
        simpleVTKWriter layerVTK
        (
            pp_.localFaces(),
            sPts
        );

        layerVTK.addPointData("pNormals", pointNormals);
        layerVTK.write("pNormSmoothed.vtk");
    }

    return pointNormals;
}

Foam::List<Foam::labelList>
Foam::edgeClassification::uniqueEdgePatches()
{
   const polyBoundaryMesh& patches = mesh_.boundaryMesh();
   List<labelList> edgePatches(pp_.nEdges(),labelList(0));

   forAll(pp_.edges(), edgei)
   {
      const labelList& edgeFaces = pp_.edgeFaces()[edgei];
      forAll(edgeFaces, eFI)
      {
         label facei = pp_.addressing()[edgeFaces[eFI]];
         label patchi = patches.whichPatch(facei);
         if (patchi != -1 && !patches[patchi].coupled())
         {
            labelList& ePatches = edgePatches[edgei];
            if (findIndex(ePatches, patchi) == -1)
            {
                label sz = edgePatches[edgei].size();
                ePatches.setSize(sz+1);
                ePatches[sz] = patchi;
            }
         }
      }
   }

   syncTools::syncEdgeList
   (
       mesh_,
       ppMeshEdges_,
       edgePatches,
       uniqueEqOp(),
       labelList()       // initial value
   );

   return edgePatches;

}

Foam::List<Foam::Tuple2<Foam::edgeClassification::edgeType,Foam::scalar>>
Foam::edgeClassification::setEdgeType
(
    const scalar convexFeatureAngle,
    const scalar concaveFeatureAngle,
    const scalar baffleFeatureAngle
)
{
    scalar featureAngle = max(convexFeatureAngle,concaveFeatureAngle);

    label nPPEdges = pp_.edges().size();
    List<pointField> edgeNormals(nPPEdges,pointField(0));

    //Check for a valid flipMap
    bool flipMapValid(flipMap_.valid() ? true : false);

    forAll(pp_.edges(), edgei)
    {
        const labelList& eFaces = pp_.edgeFaces()[edgei];
        pointField bVecs(eFaces.size());
        label nPts = 0;
        forAll(eFaces, eFI)
        {
            label lfacei = eFaces[eFI];
            label facei = pp_.addressing()[lfacei];
            bool flip(flipMapValid ? flipMap_()[lfacei] : false);
            point fA = mesh_.faces()[facei].areaNormal(pts_);

            vector fN
            (
                flip ? -fA : fA
            );
            fN /= (mag(fN) + VSMALL);
            bVecs[nPts++] = fN;
        }
        bVecs.setSize(nPts);
        edgeNormals[edgei]  = std::move(bVecs);
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        edgeNormals,
        edgeClassification::sumEqOp(),
        pointField(0)          // null value
    );

    scalarField eDotProd(nPPEdges, GREAT);
    forAll(pp_.edges(), edgei)
    {
        const pointField& eN =  edgeNormals[edgei];
        if (eN.size() == 2)
        {
            eDotProd[edgei] = min(eDotProd[edgei],(eN[0] & eN[1]));
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        eDotProd,
        minEqOp<scalar>(),
        GREAT          // null value
    );

    labelList eT(nPPEdges, -1);
    forAll(pp_, i)
    {
        label facei = pp_.addressing()[i];
        const labelList& fEdges = pp_.faceEdges()[i];
        forAll(fEdges, fEI)
        {
            label edgei = fEdges[fEI];
            label meshedgei = ppMeshEdges_[edgei];
            scalar eAngle = eDotProd[edgei];

            if (eAngle < featureAngle)
            {
                if (eAngle < -0.9902)
                {
                    if (eAngle < baffleFeatureAngle)
                    {
                        eT[edgei] = 3;
                    }
                    else if (eAngle < convexFeatureAngle)
                    {
                        eT[edgei] = 2;
                    }
                    else
                    {
                        eT[edgei] = 0;
                    }
                }
                else
                {
                    edge e = mesh_.edges()[meshedgei];
                    point fC = mesh_.faces()[facei].centre(pts_);

                    point lp0 = pts_[e[0]];
                    point lp1 = pts_[e[1]];
                    vector offset = (lp1-lp0);
                    lp0 -= offset;
                    lp1 += offset;
                    pointHit lHit = linePointRef(lp0, lp1).nearestDist(fC);

                    vector eCtofC = vector::zero;
                    if (lHit.hit())
                    {
                        eCtofC = lHit.hitPoint() - fC ;
                        eCtofC /= (mag(eCtofC) + VSMALL);
                    }
                    else
                    {
                        point eC = e.centre(pts_);
                        eCtofC = eC - fC;
                        eCtofC /= (mag(eCtofC) + VSMALL);
                    }
                    const pointField& eN =  edgeNormals[edgei];
                    scalar maxDotProd = -GREAT;
                    scalar minDotProd = GREAT;
                    forAll(eN, eNI)
                    {
                        scalar dProd = (eN[eNI] & eCtofC);

                        if (mag(dProd) > maxDotProd)
                        {
                            maxDotProd = mag(dProd);
                            minDotProd = dProd;
                        }
                    }

                    if (minDotProd > 0)
                    {
                        if (eAngle < concaveFeatureAngle)
                        {
                            eT[edgei] = 1;
                        }
                        else
                        {
                            eT[edgei] = 0;
                        }
                    }
                    else
                    {
                        if (eAngle < baffleFeatureAngle)
                        {
                            eT[edgei] = 3;
                        }
                        else if (eAngle < convexFeatureAngle)
                        {
                            eT[edgei] = 2;
                        }
                        else
                        {
                            eT[edgei] = 0;
                        }
                    }
                }
            }
            else
            {
                eT[edgei] = 0;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        eT,
        maxEqOp<label>(),
        label(-1)
    );

    List<Tuple2<edgeType,scalar>> eType
    (
        pp_.edges().size(),
        Tuple2<edgeType,scalar>(MANIFOLD,-GREAT)
    );

    forAll(ppMeshEdges_, edgei)
    {
        label nNbrs =  edgeNormals[edgei].size();
        if (nNbrs != 2)
        {
           if (nNbrs == 1)
           {
               eType[edgei] = Tuple2<edgeType,scalar>(BOUNDARY,-GREAT);
           }
           else
           {
               eType[edgei] = Tuple2<edgeType,scalar>(NONMANIFOLD,-GREAT);
           }
        }
        else
        {
            scalar eAngle = eDotProd[edgei];
            if (eT[edgei] == 0)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(MANIFOLD,eAngle);
            }
            else if (eT[edgei] == 1)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(CONCAVE,eAngle);
            }
            else if (eT[edgei] == 2)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(CONVEX,eAngle);
            }
            else if (eT[edgei] == 3)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(BAFFLE,eAngle);
            }
        }
    }

    return eType;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeClassification::edgeClassification
(
    const polyMesh& mesh,
    const pointField& pts,
    const indirectPrimitivePatch& pp,
    const labelList& ppMeshEdges,
    const scalar convexFeatureAngle,
    const scalar concaveFeatureAngle,
    const scalar baffleFeatureAngle,
    const autoPtr<boolList> flipMap
)
:
    mesh_(mesh),
    pts_(pts),
    pp_(pp),
    ppMeshEdges_(ppMeshEdges),
    flipMap_(flipMap),
    eType_
    (
        setEdgeType
        (
            convexFeatureAngle,
            concaveFeatureAngle,
            baffleFeatureAngle
        )
    )
{}

// ************************************************************************* //
