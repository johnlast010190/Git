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

#include "deformingBody/fvMeshGIBChangersDeformingBody.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "GIBTools/GIBMapping/GIBMapping.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(deformingBody, 0);
    addToRunTimeSelectionTable(fvMeshGIBChanger, deformingBody, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::deformingBody::initialization()
{
    writeInterface_ =
        dynamicMeshCoeffs_.lookupOrDefault<Switch>("writeInterface", false);

    if (dynamicMeshCoeffs_.found("motionFunctions"))
    {
        const dictionary& mfDict =
            dynamicMeshCoeffs_.subDict("motionFunctions");
        DBMFPtr_ = motionFunction::New(mesh(), mfDict, mesh().time());
    }

    mesh().faceZones().instance() = mesh().time().timeName();
    mesh().cellZones().instance() = mesh().time().timeName();

    postPro();

    Info<< dynamicMeshCoeffs_ << endl;
}


void Foam::fvMeshGIBChangers::deformingBody::makeRecAllPoints0() const
{
    if (recAllPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recAllPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recAllPoints0 = *recAllPoints0Ptr_;
    tmp<vectorField> tOldBoundary = oldBoundaryLocation();
    recAllPoints0 = tOldBoundary.ref();

    syncPoints(recAllPoints0);
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::deformingBody::movePolyPatch
(
    primitivePatch& pp,
    const vectorField& interSpeed
)
{
    PrimitivePatchInterpolation<primitivePatch> pInterC(pp);

    scalarField interSpeedmag(mag(interSpeed));
    tmp<vectorField> pinterSpeedtmp
    (
        pInterC.faceToPointInterpolate(interSpeed)
    );
    tmp<scalarField> interSpeedmagtmp
    (
        pInterC.faceToPointInterpolate(interSpeedmag)
    );

    vectorField& pinterSpeed = pinterSpeedtmp.ref();
    pinterSpeed /= (mag(pinterSpeed) + SMALL);
    scalarField& pinterSpeedmag = interSpeedmagtmp.ref();

    vectorField pnf(pinterSpeed*pinterSpeedmag);

    labelList dispCount(pp.nPoints(), 1);

    syncTools::syncPointList
    (
        mesh(),
        pp.meshPoints(),
        pnf,
        plusEqOp<point>(),
        vector::zero,
        distributionMap::transform()
    );
    syncTools::syncPointList
    (
        mesh(),
        pp.meshPoints(),
        dispCount,
        plusEqOp<label>(),
        label(0),
        distributionMap::transform()
    );

    forAll(pnf, pI)
    {
        pnf[pI] = pnf[pI]/dispCount[pI];
    }

    tmp<vectorField> newPpPointst(new vectorField(mesh().points()));
    pointField& newPpPoints = newPpPointst.ref();

    const labelList& patchPoints = pp.meshPoints();
    forAll(patchPoints, pI)
    {
        const label& ppI = patchPoints[pI];
        newPpPoints[ppI] += mesh().time().deltaTValue()*pnf[pI];
    }

    return newPpPointst;
}


void Foam::fvMeshGIBChangers::deformingBody::writeInterface
(
    const indirectPolyPatch& gibPolyPatch
)
{
    if (writeInterface_ && mesh().time().outputTime())
    {
        simpleVTKWriter
        (
            gibPolyPatch.localFaces(),
            gibPolyPatch.localPoints()
        ).write
        (
            postProFolder_
          + "/"
          + word("vis")
          + mesh().time().timeName()
          + ".vtk"
        );
    }
}


void Foam::fvMeshGIBChangers::deformingBody::postPro()
{
    if (Pstream::master() || !Pstream::parRun())
    {
        postProFolder_ =
        (
            Pstream::parRun()
          ? mesh().time().path()/".."
          : mesh().time().path()
        );
        postProFolder_ += ("/postProcessing/GIBInterface");

        mkDir(postProFolder_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::deformingBody::deformingBody(fvMesh& mesh)
:
    fvMeshGIBChangersBase(mesh, typeName),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    initialization();
}


Foam::fvMeshGIBChangers::deformingBody::deformingBody
(
    fvMesh& mesh,
    const word& typeN
)
:
    fvMeshGIBChangersBase(mesh, typeN),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    initialization();
}


Foam::fvMeshGIBChangers::deformingBody::deformingBody
(
    fvMesh& mesh,
    const dictionary dict
)
:
    fvMeshGIBChangersBase(mesh, dict),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    initialization();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::deformingBody::~deformingBody()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::deformingBody::updateInit(const word& zoneName)
{
    IOobject triIO
    (
        triName_,
        mesh().time().constant(),
        "triSurface",
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    const fileName gibSurf = mesh().time().constant()/"triSurface"/triName_;

    if (isFile(gibSurf))
    {
        ibMeshPtr_ = new triSurfaceMesh(triIO);
    }
    else
    {
        FatalErrorInFunction
            << "Could not find file " << gibSurf << nl << exit(FatalError);
    }

    fvMeshGIBChangersBase::updateInit(zoneName);
}


bool Foam::fvMeshGIBChangers::deformingBody::update()
{
    storeOldTimes();
    oldPoints_ = mesh().points();
    prevPoints_ = mesh().points();

    DBMFPtr_->update();

    const fvPatch& gibPatch(mesh().boundary()[masterGIB_]);

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    faceList faces(gibPolyPatch);
    const boolList& flipMap = gibPolyPatch.fm();

    forAll(faces, fI)
    {
        if (flipMap[fI])
        {
            faces[fI].flip();
        }
    }

    pointField pointsF(mesh().points());
    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    pointsF = movePolyPatch(pp, DBMFPtr_->interfaceVelocity(masterGIB_)());

    // Important to update pp (cleaning geometry)
    pp.clearGeom();

    GIBMapping mapCl(*this, pp);

    //- Careful: if this is destroyed, ibMeshPtr_ reference is wrong
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
                IOobject::NO_WRITE,
                false
            ),
            mapCl.triS()
        );

    writeInterface(gibPolyPatch);

    clearOutGIBData();

    doUpdate(mapCl, false);

    return true;
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::deformingBody::velocityCorrect
(
    const vectorField& pc
) const
{
    return DBMFPtr_->boundaryVelocity(masterGIB_);
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::deformingBody::oldBoundaryLocation() const
{
    return DBMFPtr_->interfacePointsVelocity(masterGIB_);
}


// ************************************************************************* //
