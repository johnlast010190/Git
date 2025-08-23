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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "adjoint/fvMeshGIBChangersAdjoint.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "GIBTools/GIBAreaSmoothing/GIBAreaSmoothing.H"
#include "GIBTools/GIBMapping/GIBMapping.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(adjoint, 0);
    addToRunTimeSelectionTable(fvMeshGIBChanger, adjoint, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::adjoint::courantScaling(scalarField& pG)
{
    const labelList& patchPoints =
        mesh().boundary()[masterGIB_].patch().meshPoints();
    const labelListList& pEdges = mesh().pointEdges();

    scalarField adgLp = scalarField(patchPoints.size(), 0);

    forAll(patchPoints, pI)
    {
        label pointI = patchPoints[pI];
        const labelList& pEdgesI = pEdges[pointI];
        forAll(pEdgesI, eI)
        {
            label edgeI = pEdgesI[eI];
            const edge& edgeII = mesh().edges()[edgeI];
            label point1 = edgeII.start();
            label point2 = edgeII.end();
            scalar edgL =
                mag(mesh().points()[point1] - mesh().points()[point2]);
            adgLp[pI] += edgL;
        }
        adgLp[pI] /= pEdgesI.size();
    }

    scalarField pCourant = scalarField(patchPoints.size(), 0);

    scalar maxCourant =
        adjProperties().subDict("topology")
            .subDict("surfaceTracking")
            .lookupOrDefault<scalar>("maxCourant", 0.5);

    autoPtr<Function1<scalar>> meanCourantPtr
        (
            Function1<scalar>::New
            (
                "meanCourant",
                adjProperties().subDict("topology").subDict("surfaceTracking")
            )
        );

    scalar meanCourant(meanCourantPtr->value(mesh().time().timeOutputValue()));

    for (int i = 1; i < 20; i++)
    {
        pCourant = mag(pG)/adgLp;
        scalar meanC = gSum(mag(pG))/gSum(adgLp);

        if (meanC != 0)
        {
            scalar scale = meanCourant/meanC;
            pG *= scale;

            forAll(pG, pI)
            {
                scalar maCI = mag(pG[pI])/adgLp[pI];
                if (maCI > maxCourant)
                {
                    pG[pI] = maxCourant*sign(pG[pI])*adgLp[pI];
                }
            }
            pCourant = mag(pG)/adgLp;
        }
    }

    Info<< endl;
    Info<< "Courant max: "  << gMax(pCourant) << endl;
    Info<< "Courant mean: " << gSum(mag(pG))/gSum(adgLp) << endl;
}


Foam::scalarField Foam::fvMeshGIBChangers::adjoint::calcCurvature()
{
    scalarField curv = mesh().boundary()[masterGIB_].patch().pointCurvature();
    scalarField curvTmp = curv;

    const labelListList& pEdges =
        mesh().boundary()[masterGIB_].patch().pointEdges();
    const edgeList& edges = mesh().boundary()[masterGIB_].patch().edges();

    forAll(pEdges, pI)
    {
        const labelList& pEdge = pEdges[pI];
        forAll(pEdge, eI)
        {
            label edgeI = pEdge[eI];
            const edge& edgeII = edges[edgeI];
            label point1 = edgeII.start();
            label point2 = edgeII.end();

            if (point1 != pI)
            {
                curvTmp[pI] += curv[point1];
            }
            else
            {
                curvTmp[pI] += curv[point2];
            }
        }
        curvTmp[pI] /= pEdges[pI].size();
    }

    curv = curvTmp;

    return curv;
}


void Foam::fvMeshGIBChangers::adjoint::filterBoundaryPoints
(
    scalarField& sensN,
    const indirectPolyPatch& gibPolyPatch
)
{
    const pointField& basePoints = *basePoints_;
    const label& zoneId = gibFaceZone();
    const faceZone& flist = mesh().faceZones()[zoneId];

    forAll(flist, fI)
    {
        label flistI = flist[fI];
        if (flistI >= mesh().nInternalFaces())
        {
            forAll(mesh().faces()[flistI], pI)
            {
                label gpointI = mesh().faces()[flistI][pI];
                label pointI = gibPolyPatch.whichPoint(gpointI);
                if (pointI != -1)
                {
                    if
                    (
                        mesh().points()[gpointI] == basePoints[gpointI]
                     && sensN[pointI] > 0
                    )
                    {
                        sensN[pointI] = 0;
                    }
                }
            }
        }
    }
}


Foam::boolList Foam::fvMeshGIBChangers::adjoint::findConstraintFaces()
{
    const pointField& basePoints = *basePoints_;
    const volScalarField& G = mesh().lookupObject<volScalarField>("G");
    const polyPatch& poly = mesh().boundary()[masterGIB_].patch();
    boolList conFaces(poly.size(), false);

    if (isA<indirectPolyPatch>(poly))
    {
        const indirectPolyPatch& inPoly =
            refCast<const indirectPolyPatch>(poly);

        const labelList& addr = inPoly.fAddr();
        forAll(addr, fI)
        {
            label addrI = addr[fI];
            if (addrI >= mesh().nInternalFaces())
            {
                bool moved = false;
                forAll(mesh().faces()[addrI], pI)
                {
                    label gpointI = mesh().faces()[addrI][pI];
                    if
                    (
                        (mesh().points()[gpointI] != basePoints[gpointI])
                     && !moved
                    )
                    {
                        moved = true;
                    }
                }
                if ((!moved) && (G.boundaryField()[masterGIB_][fI] > 0))
                {
                    conFaces[fI] = true;
                }
            }
        }
    }

    return conFaces;
}


void Foam::fvMeshGIBChangers::adjoint::fixConstraintPatches(pointField& snapP)
{
    const pointField& basePoints = *basePoints_;
    forAll(mesh().boundary(), pI)
    {
        if
        (
            mesh().boundary()[pI].type() == "inlet"
         || mesh().boundary()[pI].type() == "outlet"
         || mesh().boundary()[pI].patch().physicalType() == "inlet"
         || mesh().boundary()[pI].patch().physicalType() == "outlet"
        )
        {
            const labelList& pP = mesh().boundary()[pI].patch().meshPoints();
            forAll(pP, pI)
            {
                snapP[pP[pI]] = basePoints[pP[pI]];
            }
        }
    }
}


void Foam::fvMeshGIBChangers::adjoint::syncProcBoundaryPoints
(
    pointField& newPoints,
    const pointField& newPointsCpy
)
{
    const labelList& patchPoints =
        mesh().boundary()[masterGIB_].patch().meshPoints();

    pointField relDis(newPointsCpy - newPoints);

    syncTools::syncPointList
    (
        mesh(),
        patchPoints,
        relDis,
        plusEqOp<vector>(),
        vector::zero
    );

    newPoints += relDis;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::fvMeshGIBChangers::adjoint::computeNewPoints
(
    primitivePatch& pp,
    const scalarField& interSpeed
)
{
    PrimitivePatchInterpolation<primitivePatch> pInterC(pp);

    tmp<Field<scalar>> pinterSpeedtmp
    (
        new Field<scalar>(pInterC.faceToPointInterpolate(interSpeed))
    );

    scalarField& pinterSpeed = pinterSpeedtmp.ref();

    const fvPatch& gibPatch(mesh().boundary()[masterId()]);
    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    if (false)
    {
        simpleVTKWriter a1(pp.localFaces(), pp.localPoints());

        a1.addPointData("pinterSpeed", pinterSpeed);
    }

    autoPtr<Function1<scalar>> curvWeightPtr
    (
        Function1<scalar>::New
        (
            "curvatureWeight",
            adjProperties().subDict("topology").subDict("surfaceTracking")
        )
    );

    scalar curvWeight(curvWeightPtr->value(mesh().time().timeOutputValue()));

    scalarField curv(curvWeight*calcCurvature());

    tmp<Field<vector>> ppnftmp
    (
        new Field<vector>
        (
            pInterC.faceToPointInterpolate(mesh().boundary()[masterGIB_].nf())
        )
    );
    const vectorField& np = ppnftmp();

    scalarField pG(pinterSpeed - curv);

    filterBoundaryPoints(pG, gibPolyPatch);

    courantScaling(pG);

    vectorField pnf(pG*np);

    // Syncing the processor boundary points
    pointField pf1 = pointField(mesh().points().size(), vector::zero);
    pointField pf2 = pointField(mesh().points().size(), vector::zero);
    const labelList& patchPoints = gibPolyPatch.meshPoints();
    forAll(patchPoints, pI)
    {
        label ppI = patchPoints[pI];
        pf1[ppI] = pnf[pI];
        pf2[ppI] = pnf[pI];
    }

    applyTwoDPlanesCorrection(pf1, pf2);

    syncTools::syncPointPositions
    (
        mesh(),
        pf1,
        maxEqOp<point>(),              // combine op
        point(-GREAT,-GREAT,-GREAT)    // null
    );
    syncTools::syncPointPositions
    (
        mesh(),
        pf2,
        minEqOp<point>(),           // combine op
        point(GREAT,GREAT,GREAT)    // null
    );

    forAll(patchPoints, pI)
    {
        const label& ppI = patchPoints[pI];
        pnf[pI] = (pf1[ppI]+pf2[ppI])/2;
    }

    tmp<Field<vector>> newPpPointst(new Field<vector>(pp.points()));
    pointField& newPpPoints = newPpPointst.ref();

    forAll(newPpPoints, pI)
    {
        newPpPoints[pI] += mesh().time().deltaTValue()*pnf[pI];
    }

    return newPpPointst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::adjoint::adjoint(fvMesh& mesh)
:
    deformingBody(mesh),
    adjPropertiesPtr_(nullptr)
{}


Foam::fvMeshGIBChangers::adjoint::adjoint(fvMesh& mesh, const dictionary dict)
:
    deformingBody(mesh, dict),
    adjPropertiesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::adjoint::~adjoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::adjoint::updateInit()
{
    findGIBPatches();
    clearOutGIBData();

    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = mesh().faceZones()[zoneId];
    faceZone& fZone = const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    DynamicList<label> dfl(mesh().faces().size());
    forAll(mesh().boundary(), pI)
    {
        if (isA<wallFvPatch>(mesh().boundary()[pI]))
        {
            label startPI = mesh().boundary()[pI].start();
            forAll(mesh().boundary()[pI], pfI)
            {
                dfl.append(startPI+pfI);
            }
        }
    }
    dfl.shrink();

    fZone.resize(dfl.size());
    fZone.resetAddressing(dfl, boolList(dfl.size(), false));

    mesh().updateGIB();
}


bool Foam::fvMeshGIBChangers::adjoint::update()
{
    storeOldTimes();
    oldPoints_ = mesh().points();

    //----------------------------------------------------------//
    DBMFPtr_->update();

    const fvPatch& gibPatch(mesh().boundary()[masterId()]);
    faceList faces = preparePatch(gibPatch);

    pointField pointsF = gibPatch.patch().localPoints();

    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    pointField oldp = pointsF;

    GIBAreaSmoothing saSmoother
    (
        mesh(),
        adjProperties().subDict("topology").subDict("surfaceTracking"),
        masterGIB_,
        findConstraintFaces()
    );

    saSmoother.update();

    tmp<vectorField> newPpPointst =
        computeNewPoints(pp, saSmoother.smoothSens()());

    pointsF = newPpPointst();

    nearBoundaryIntersectionsChecking(pointsF);

    checkConcaveBoundaryPatchPoints(pointsF);

    GIBMapping mapCl = GIBMapping(*this, pp, false);

    deleteDemandDrivenData(ibMeshPtr_);
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                mesh().time().constant(),
                "triSurface",
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mapCl.triS()
        );

    clearOutGIBData();
    //----------------------------------------------------------//

    const label& zoneId = gibFaceZone();
    const faceZone& cfZone = mesh().faceZones()[zoneId];
    faceZone& fZone = const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    tmp<pointField> snapPt = findSnappedPoints();
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList&  fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    correctConstraintPatches(snapP);

    fixConstraintPatches(snapP);
    correctBoundaryPointsOnBaseMesh(snapP);

    syncPoints(snapP);

    mesh().moveGIBPoints(snapP);
    mesh().updateGIB();

    faceCellsVisDebug();

    mapCl.mapBcs();

    popShrinkFields();
    correctV0();
    popGrowFields();
    resetMeshFluxes();
    correctBCs();

    return true;
}


const Foam::IOdictionary& Foam::fvMeshGIBChangers::adjoint::adjProperties()
{
    if (!adjPropertiesPtr_)
    {
        adjPropertiesPtr_ =
        &(
            mesh().time().lookupObject<IOdictionary>("adjointProperties")
        );
    }

    return *adjPropertiesPtr_;
}


// ************************************************************************* //
