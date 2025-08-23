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

#include "solidBody/fvMeshGIBChangersSolidBody.H"
#include "GIBTools/GIBMapping/GIBMapping.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(solidBody, 0);
    addToRunTimeSelectionTable(fvMeshGIBChanger, solidBody, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::solidBody::initialise()
{
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                mesh().time().constant(),
                "triSurface",
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    initStlPointsPtr_ = new pointField(ibMeshPtr_->points());

    if (dynamicMeshCoeffs_.found("motionFunctions"))
    {
        const dictionary& motionFDict =
            dynamicMeshCoeffs_.subDict("motionFunctions");
        isSolidBody_ = motionFDict.found("solidBodyMotionFunction");

        if (isSolidBody_)
        {
            SBMFPtr_ = solidBodyMotionFunction::New(motionFDict, mesh().time());

            DeprecationWarningInFunction
            (
                solidBodyMotionFunction::typeName,
                "mesh motion type",
                40200,
                "Please replace it by using referenceFrame."
            );
        }
        else
        {
            coorFramePtr_ = coordinateFrame::lookupNew(mesh(), motionFDict);
            coorFramePtr_->resetDynamic(true);
            // Future check for dynamic motion
            if (!coorFramePtr_->anyDynamic())
            {
                FatalErrorInFunction
                    << "coordinateFrame " << coorFramePtr_->name()
                    << " is linked to dynamic mesh but it is not dynamic."
                    << abort(FatalError);
            }
        }
    }

    if (dynamicMeshCoeffs_.found("geometryTransformationFunctions"))
    {
        const dictionary& geoTDict =
            dynamicMeshCoeffs_.subDict("geometryTransformationFunctions");

        geoTrans_.setSize(geoTDict.size());

        label fID = 0;
        forAllConstIter(dictionary, geoTDict, iter)
        {
            if (iter().isDict())
            {
                const dictionary& subDict = iter().dict();

                geoTrans_.set
                (
                    fID,
                    geometryTransformation::New(subDict)
                );
                fID++;
            }
        }
        geoTrans_.setSize(fID);
    }

    if (dynamicMeshCoeffs_.found("solverMotionFunctions"))
    {
        const dictionary& solverFDict =
            dynamicMeshCoeffs_.subDict("solverMotionFunctions");
        motionPtr_ = motionSolver::New(mesh(), solverFDict);
    }

    // Move in construction to store the old points,
    // otherwise the old points will be the original points of stl,
    // In construction it is using fvMeshGIBChangersSolidBody::moveSurface,
    // however later it uses virtual dispatch.
    tmp<pointField> newSurfPointst = moveSurface(false);
    pointField& newSurfPoints = newSurfPointst.ref();

    ibMeshPtr_->movePoints(newSurfPoints);

    mesh().faceZones().instance() = mesh().time().timeName();
    mesh().cellZones().instance() = mesh().time().timeName();

    calculateHitIndex();

    Info<< dynamicMeshCoeffs_ << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::fvMeshGIBChangers::solidBody::movePolyPatch
(
    primitivePatch& pp,
    const vectorField& interSpeed
)
{
    tmp<vectorField> newPpPointst(new vectorField(pp.points()));
    pointField& newPpPoints = newPpPointst.ref();
    forAll(newPpPoints, pI)
    {
        newPpPoints[pI] += mesh().time().deltaTValue()*interSpeed[pI];
    }

    return newPpPointst;
}


Foam::tmp<Foam::pointField>
Foam::fvMeshGIBChangers::solidBody::moveSurface(bool updateSolver)
{
    tmp<pointField> stlPointst(new pointField(ibMeshPtr_->points()));
    pointField& stlPoints = stlPointst.ref();
    const vectorField& initStlPoints = *initStlPointsPtr_;

    if (SBMFPtr_.valid() || coorFramePtr_)
    {
        septernion transform =
            coorFramePtr_
          ? coorFramePtr_->transformation()
          : SBMFPtr_().transformation();

        if (Pstream::parRun())
        {
            List<septernion> transformList(Pstream::nProcs());
            transformList[Pstream::myProcNo()] = transform;

            // Distribute transformation
            Pstream::allGatherList(transformList);

            // Do transformation based on master processor
            transform = transformList[Pstream::masterNo()];
        }

        const bool isIncremental =
            (coorFramePtr_ && coorFramePtr_->isIncrementalMotion())
         || (SBMFPtr_.valid() && SBMFPtr_().isIncrementalMotion());

        if (isIncremental)
        {
            // Puting the initialisation of old points here should deal even
            // with the delayed construction of some frames where
            // the incremental flag could be set later.
            if (!geomOldPoints_.valid())
            {
                geomOldPoints_.reset(new pointField(initStlPoints));
            }

            stlPoints = transformPoints(transform, geomOldPoints_());
            geomOldPoints_() = stlPoints;
        }
        else
        {
            stlPoints = transformPoints(transform, initStlPoints);
        }
    }

    pointField oldPoints = stlPoints;

    forAll(geoTrans_, i)
    {
        stlPoints = geoTrans_[i].transformPoints(oldPoints)();
    }

    if (motionPtr_.valid())
    {
        if (updateSolver)
        {
            stlPoints = motionPtr_->newPoints();
        }
        else
        {
            stlPoints = motionPtr_->curPoints();
        }
    }

    return stlPointst;
}


Foam::tmp<Foam::pointField>
Foam::fvMeshGIBChangers::solidBody::recNewPointLocation
(
    const pointField& pf0,
    const labelList& addr
)
{
    tmp<pointField> pft(new pointField(pf0.size(), Zero));
    pointField& pf = pft.ref();

    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    const labelList& hitInd = hitIndex();

    // Size checks for debugging
    if (addr.size() != pf0.size())
    {
        FatalErrorInFunction
            << "indexes of the points and the pointField "
            << "are not the same size"
            << abort(FatalError);
    }

    if (mesh().points().size() != hitInd.size())
    {
        if (hitInd.size() != 0)
        {
            FatalErrorInFunction
                << "Inconsistent sizes between hitIndices and points"
                << abort(FatalError);
        }
        else
        {
            Info<< "hitIndex is not set." << nl
                << "triSurface constructed from a faceZone." << nl
                << "Hit index is created on the fly."
                << endl;
        }
    }

    // We want to express the point based on the coordinates of the
    // points of the hit triangle of the stl.
    // p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a 3x3 system
    forAll(addr, pI)
    {
        label addrI = addr[pI];
        label hitIndexI = hitInd[addrI];
        const triFace& triList = triSM().triSurface::operator[](hitIndexI);

        vector p0 = pf0[pI];

        vector p_1 = surP[triList[0]];
        vector p_2 = surP[triList[1]];
        vector p_3 = surP[triList[2]];

        vector w = triList.findTriangleWeights(p0, surP0);

        pf[pI] = w[0]*p_1 + w[1]*p_2 + w[2]*p_3;
    }

    return pft;
}


void Foam::fvMeshGIBChangers::solidBody::computeOldPositionsInUnsnappedCells
(
    pointField& recAllPoints0,
    const labelList& fscp
) const
{
    pointField currentPoints(fscp.size(), Zero);
    forAll(fscp, pI)
    {
        label gpI = fscp[pI];
        currentPoints[pI] = mesh().points()[gpI];
    }

    List<pointIndexHit> nearest;
    triSM().findNearest
    (
        currentPoints,
        scalarField(currentPoints.size(), sqr(GREAT)),
        nearest
    );
    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    forAll(currentPoints, pI)
    {
        label hitIndexI = nearest[pI].index();
        const triFace& triList = triSM().triSurface::operator[](hitIndexI);
        vector disp = surP0[triList[0]]-surP[triList[0]];

        recAllPoints0[fscp[pI]] = currentPoints[pI]+disp;
    }
}


void Foam::fvMeshGIBChangers::solidBody::calculateHitIndex()
{
    labelIOList& hitIndex = *hitIndexPtr_;

    if (mesh().points().size() != hitIndex.size())
    {
        if (hitIndex.size() != 0)
        {
            FatalErrorInFunction
                << "Inconsistent sizes between hitIndeces and points"
                << abort(FatalError);
        }
        else
        {
            Info<< "hitIndex is not set."
                << "triSurface constructed from a faceZone."
                << "hit Index is created on the fly."
                << endl;

            if (hitIndex.size()==0)
            {
                hitIndex.resize(mesh().points().size(), -1);

                List<pointIndexHit> nearest;
                triSM().findNearest
                (
                    mesh().points(),
                    scalarField(mesh().points().size(), sqr(GREAT)),
                    nearest
                );
                forAll(nearest, pI)
                {
                    const pointIndexHit& nearI = nearest[pI];
                    hitIndex[pI] = nearI.index();
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::solidBody::solidBody(fvMesh& mesh)
:
    fvMeshGIBChangersBase(mesh, typeName),
    initStlPointsPtr_(nullptr),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    initialise();
}


Foam::fvMeshGIBChangers::solidBody::solidBody(fvMesh& mesh, const word& typeN)
:
    fvMeshGIBChangersBase(mesh, typeN),
    initStlPointsPtr_(nullptr),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                mesh.time().constant(),
                "triSurface",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    initStlPointsPtr_ = new pointField(ibMeshPtr_->points());

    mesh.faceZones().instance() = mesh.time().timeName();
    mesh.cellZones().instance() = mesh.time().timeName();

    Info<< dynamicMeshCoeffs_ << endl;
}


Foam::fvMeshGIBChangers::solidBody::solidBody
(
    fvMesh& mesh,
    const dictionary dict
)
:
    fvMeshGIBChangersBase(mesh, dict),
    initStlPointsPtr_(nullptr),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr),
    geoTrans_(0),
    motionPtr_()
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::solidBody::~solidBody()
{
    delete initStlPointsPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::solidBody::updateInit(const word& zoneName)
{
    tmp<pointField> newSurfPoints = moveSurface(false);

    ibMeshPtr_->movePoints(newSurfPoints());

    fvMeshGIBChangersBase::updateInit(zoneName);
}


bool Foam::fvMeshGIBChangers::solidBody::update()
{
    storeOldTimes();
    oldPoints_ = mesh().points();
    prevPoints_ = mesh().points();

    if (coorFramePtr_)
    {
        coorFramePtr_->updateState();
    }

    tmp<pointField> newSurfPoints = moveSurface(true);

    ibMeshPtr_->movePoints(newSurfPoints());

    const fvPatch& gibPatch(mesh().boundary()[masterId()]);
    faceList faces = preparePatch(gibPatch);

    pointField pointsF = gibPatch.patch().localPoints();

    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    tmp<vectorField> newPpPointst =
        recNewPointLocation(pointsF, gibPatch.patch().meshPoints());

    pointsF = newPpPointst();

    pp.clearGeom();

    GIBMapping mapCl = GIBMapping(*this, pp);

    clearOutGIBData();

    doUpdate(mapCl, true);

    return true;
}


Foam::tmp<Foam::vectorField> Foam::fvMeshGIBChangers::solidBody::velocityCorrect
(
    const vectorField& pc
) const
{
    FatalErrorInFunction
        << "Boundary condition of velocity should be computed based on the "
        << "old positions of the points."
        << abort(FatalError);

    return tmp<vectorField>(new vectorField(pc.size(), Zero));
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::solidBody::oldBoundaryLocation() const
{
    const boolList& intP = interP();
    const vectorField& cP = mesh().points();

    tmp<vectorField> cP0t(new vectorField(cP.size(), Zero));
    vectorField& cP0 = cP0t.ref();

    cP0 = velocityCorrect(cP)();

    boolList oldnewintPoints = boolList(mesh().points().size(), false);
    boolList popUpPoints = boolList(mesh().points().size(), false);
    boolList marknewFaces = boolList(mesh().faces().size(), false);
    boolList markoldFaces = boolList(mesh().faces().size(), false);

    const labelList& fZone0 = fl0();
    const labelList& fZone = fl();

    forAll(fZone0, fI)
    {
        markoldFaces[fZone0[fI]] = true;
    }
    forAll(fZone, fI)
    {
        marknewFaces[fZone[fI]] = true;
    }
    forAll(mesh().faces(), fI)
    {
        if (marknewFaces[fI] && markoldFaces[fI])
        {
            forAll(mesh().faces()[fI], pI)
            {
                label pointI = mesh().faces()[fI][pI];
                oldnewintPoints[pointI] = true;
            }
        }
    }

    const boolList& popUpC = popUpCells();

    forAll(popUpC, cI)
    {
        if (popUpC[cI])
        {
            const labelList& cP = mesh().cellPoints()[cI];
            forAll(cP, pI)
            {
                popUpPoints[cP[pI]] = true;
            }
        }
    }

    boolList boundaryP(mesh().points().size(), false);
    forAll(mesh().boundary(), pI)
    {
        const polyPatch& pp = mesh().boundary()[pI].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if
            (
                !isA<emptyFvPatch>(mesh().boundary()[pI])
             && !isA<wedgeFvPatch>(mesh().boundary()[pI])
            )
            {
                if (pp.size())
                {
                    const labelList& pPoints = pp.meshPoints();
                    forAll(pPoints, pI)
                    {
                        label gpI = pPoints[pI];
                        boundaryP[gpI] = true;
                    }
                }
            }
        }
    }

    forAll(intP, pI)
    {
        if (intP[pI])
        {
            if (popUpPoints[pI])
            {
                // Previous location, if point was and is at the interface
                cP0[pI] = prevPoints_[pI];
            }
            else
            {
                // If point was not and it is at the interface,
                // calculate old location based on the interface velocity.
                cP0[pI] = cP[pI] - cP0[pI]*mesh().time().deltaTValue();
            }
        }
        else
        {
            // For all the non-interface points
            cP0[pI] = cP[pI];
        }
    }

    return cP0t;
}


// ************************************************************************* //
