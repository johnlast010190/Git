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

#include "cfdTools/general/include/fvCFD.H"
#include "meshes/polyMesh/polyPatches/derived/inlet/inletPolyPatch.H"
#include "meshes/polyMesh/polyPatches/derived/outlet/outletPolyPatch.H"
#include "meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "regionSplit/regionSplit.H"
#include "GIBTools/GIBMapping/GIBMapping.H"
#include "MeshedSurface/MeshedSurfaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshGIBChangersBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangersBase::initialize()
{
    modifyRegion0Patch();

    findGIBPatches();
    IOobject ioBasPoints
    (
        "basePoints",
        mesh().time().findInstance
        (
            mesh().meshDir(),
            "basePoints",
            IOobject::READ_IF_PRESENT
        ),
        mesh().meshSubDir,
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (exists(ioBasPoints.objectPath()))
    {
        basePoints_ = new pointIOField(ioBasPoints);
    }
    else
    {
        basePoints_ =
            new pointIOField
            (
                IOobject
                (
                    "basePoints",
                    mesh().polyMesh::instance(),
                    mesh().meshSubDir,
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );
        pointField& basePoints = *basePoints_;
        basePoints = mesh().points();
        basePoints_->write();
    }

    IOobject ioHitIndex
    (
        "hitIndex",
        mesh().time().findInstance
        (
            mesh().meshDir(),
            "hitIndex",
            IOobject::READ_IF_PRESENT
        ),
        mesh().meshSubDir,
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (exists(ioHitIndex.objectPath()))
    {
        hitIndexPtr_ = new labelIOList(ioHitIndex);
    }
    else
    {
        hitIndexPtr_ =
            new labelIOList
            (
                IOobject
                (
                    "hitIndex",
                    mesh().polyMesh::instance(),
                    mesh().meshSubDir,
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }

    baseSf_ = new vectorField(mesh().faces().size(), Zero);
    baseCf_ = new vectorField(mesh().faces().size(), Zero);
    vectorField& baseSf = *baseSf_;
    vectorField& baseCf = *baseCf_;

    scalarField magfAreas(mesh().faces().size());

    baseCC_ = new vectorField(mesh().cells().size(), Zero);
    vectorField& baseCC = *baseCC_;

    scalarField cellVols(scalarField(baseCC.size(), 0));
    const pointField& basePoints = *basePoints_;

    mesh().makeFaceCentresAndAreas(basePoints, baseCf, baseSf, magfAreas);
    mesh().makeCellCentresAndVols(baseCf, baseSf, baseCC, cellVols);
    makeConcavePoints();
    makeNormal2DEdges();

    mesh().faceZones().instance() = mesh().time().timeName();
    mesh().cellZones().instance() = mesh().time().timeName();
    mesh().faceZones().writeOpt() = IOobject::AUTO_WRITE;
    mesh().cellZones().writeOpt() = IOobject::AUTO_WRITE;

    const label zoneId = gibFaceZone();
    const faceZone& cfZone = mesh().faceZones()[zoneId];
    fl0Ptr_ = new labelList(cfZone);
    fm0Ptr_ = new boolList(cfZone.flipMap());

    Info<< dynamicMeshCoeffs_ << endl;
}


void Foam::fvMeshGIBChangersBase::modifyRegion0Patch()
{
    if (region0Patch_.size() == 0)
    {
        Info<< "Automatic search of the patch attached to fluid: " << endl;

        label patchi = 0;
        do
        {
            const polyPatch& pp = mesh().boundary()[patchi].patch();
            if
            (
                isA<inletPolyPatch>(pp)
             || isA<outletPolyPatch>(pp)
             || pp.physicalType() == "inlet"
             || pp.physicalType() == "outlet"
             || pp.physicalType() == "opening"
            )
            {
                region0Patch_.append(pp.name());
            }
            else if (isA<processorPolyPatch>(pp))
            {
                break;
            }
            patchi ++;
        } while (patchi < mesh().boundary().size());
    }
}


void Foam::fvMeshGIBChangersBase::makeFl() const
{
    if (flPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    DynamicList<label> dfl(mesh().faces().size());

    labelList interestFaces = checkingIntersectingFaces();
    const vectorField& baseCc = *baseCC_;
    const vectorField& baseCf = *baseCf_;

    pointField onCC(interestFaces.size(), Zero);
    pointField nbCC(interestFaces.size(), Zero);

    boolList tmpMultInterFaces(mesh().faces().size(), false);

    pointField neiCc(mesh().nFaces() - mesh().nInternalFaces(), Zero);
    pointField ownCc(mesh().nFaces() - mesh().nInternalFaces(), Zero);

    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& pp = mesh().boundary()[patchi].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            const labelUList& faceCells = pp.faceCells();
            if (pp.coupled())
            {
                forAll(pp, pFacei)
                {
                    const label facei = pp.start() + pFacei;
                    neiCc[facei - mesh().nInternalFaces()] =
                        baseCc[faceCells[pFacei]];
                    ownCc[facei - mesh().nInternalFaces()] =
                        baseCc[faceCells[pFacei]];
                }
            }
            else
            {
                forAll(pp, pFacei)
                {
                    const label facei = pp.start() + pFacei;
                    neiCc[facei - mesh().nInternalFaces()] = baseCf[facei];
                    ownCc[facei - mesh().nInternalFaces()] =
                        baseCc[faceCells[pFacei]];
                }
            }
        }
    }

    syncTools::swapBoundaryFacePositions(mesh(), neiCc);

    forAll(interestFaces, facei)
    {
        label interestFacesI = interestFaces[facei];
        if (interestFacesI < mesh().nInternalFaces())
        {
            const label own = mesh().owner()[interestFacesI];
            const label nei = mesh().neighbour()[interestFacesI];
            onCC[facei] = baseCc[own];
            nbCC[facei] = baseCc[nei];
        }
        else
        {
            const label patchi =
                mesh().boundaryMesh().whichPatch(interestFacesI);
            const polyPatch& pp = mesh().boundary()[patchi].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                onCC[facei] = ownCc[interestFacesI - mesh().nInternalFaces()];
                nbCC[facei] = neiCc[interestFacesI - mesh().nInternalFaces()];
            }
        }
    }

    extendBoundaryBaseMesh(interestFaces, nbCC);

    List<List<pointIndexHit>> nearest;

    if (false)
    {
        const word timeName(mesh().time().timeName());
        simpleVTKWriter(onCC).write("onCC_" + timeName + ".vtk");
        simpleVTKWriter(nbCC).write("nbCC_" + timeName + ".vtk");
    }

    triSM().findLineAll(onCC, nbCC, nearest);

    boolList blockedFace(mesh().nFaces(), false);

    forAll(interestFaces, facei)
    {
        label interestFacesI = interestFaces[facei];

        const List<pointIndexHit>& lnearI = nearest[facei];
        if (lnearI.size() > 0)
        {
            if (!(lnearI.size() % 2 == 0))
            {
                const pointIndexHit& nearI = lnearI[0];
                if (nearI.hit())
                {
                    dfl.append(interestFacesI);
                    blockedFace[interestFacesI] = true;
                }
            }
            else if (lnearI.size() == 2)
            {
                scalar d1 = mag(onCC[facei] - nbCC[facei]);
                scalar d2 = mag(lnearI[0].hitPoint() - lnearI[1].hitPoint());
                vector n1 = triSM().faceNormals()[lnearI[0].index()];
                vector n2 = triSM().faceNormals()[lnearI[1].index()];

                if (lnearI[0].hitPoint() == lnearI[1].hitPoint())
                {
                    dfl.append(interestFacesI);
                    blockedFace[interestFacesI] = true;
                }
                else
                {
                    // For multiple hiting
                    if (d2 < SMALL)
                    {
                        if ((n1&n2) > 0)
                        {
                            dfl.append(interestFacesI);
                            blockedFace[interestFacesI] = true;
                        }
                    }
                    else
                    {
                        const face& face1 =
                            triSM().triSurface::operator[](lnearI[0].index());
                        const face& face2 =
                            triSM().triSurface::operator[](lnearI[1].index());
                        bool tolIssue =
                            checkIfEdgeToleranceIntersecting(face1, face2);
                        if (tolIssue && d2 < (1e-10*d1))
                        {
                            dfl.append(interestFacesI);
                            blockedFace[interestFacesI] = true;
                        }
                    }
                }
                tmpMultInterFaces[interestFacesI] = true;
            }
            else
            {
                if (debugMode_)
                {
                    Pout<< "Intersection checking for face: "
                         << interestFacesI << tab
                         << "intersection info: " << endl;
                    Pout<< lnearI << endl;
                }
            }
        }
    }

    dfl.shrink();
    flPtr_ = new labelList(dfl);

    if (false)
    {
        const pointField& basePoints = *basePoints_;
        simpleVTKWriter
        (
            mesh().faces(),
            labelList(dfl),
            basePoints
        ).write
        (
            "fl0" + mesh().time().timeName() + ".vtk"
        );
    }

    multInterFacesPtr_ = new boolList(tmpMultInterFaces);
}


void Foam::fvMeshGIBChangersBase::makeFlipMap() const
{
    if (fmPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    calculateFlipMap();

    if (debugMode_)
    {
        regionVisDebug();
    }
}


void Foam::fvMeshGIBChangersBase::makeInterPoints() const
{
    if (interPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    interPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& interPoints = *interPointsPtr_;

    flipCells(interPoints, true);
}


void Foam::fvMeshGIBChangersBase::makeInterPoints0() const
{
    if (interPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    interPoints0Ptr_ = new boolList(mesh().points().size(), false);

    boolList& interPoints0 = *interPoints0Ptr_;
    const labelList& fll = fl0();
    forAll(fll, facei)
    {
        const face& f = mesh().faces()[fll[facei]];
        forAll(f, pFacei)
        {
            interPoints0[f[pFacei]] = true;
        }
    }
}


void Foam::fvMeshGIBChangersBase::makePhiB() const
{
    if (phiBPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    const pointField& cP = mesh().points();

    dimensionedScalar deltaT = mesh().time().deltaT();

    const vectorField& cP0 = recAllPoints0();

    tmp<scalarField> tsweptVols(new scalarField(mesh().faces().size()));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(mesh().faces(), facei)
    {
        sweptVols[facei] = mesh().faces()[facei].sweptVol(cP0, cP);
    }
    phiBPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "mPhi",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimVolume/dimTime
        );
    surfaceScalarField& mPhi = *phiBPtr_;

    scalarField& mPhiVol = mPhi.ref();
    surfaceScalarField::Boundary& mPhiB = mPhi.boundaryFieldRef();

    const fvPatchList& patches = mesh().boundary();

    mPhiVol = scalarField::subField(sweptVols, mesh().nInternalFaces());
    mPhiVol /= deltaT.value();

    // Needs recheck
    forAll(patches, patchi)
    {
        if (!isA<indirectPolyPatch>(patches[patchi].patch()))
        {
            mPhiB[patchi] = patches[patchi].patchSlice(sweptVols);
            mPhiB[patchi] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, facei)
    {
        const label fII = fzNew[facei];
        if (fII < mesh().nInternalFaces())
        {
            mPhiB[masterGIB_][facei] = mPhi[fII];
            mPhiB[slaveGIB_][facei] = -mPhi[fII];
            if (fmNew[facei])
            {
                mPhiB[masterGIB_][facei] *= -1;
                mPhiB[slaveGIB_][facei] *= -1;
            }
            fiBm += 1;
            fiBs += 1;
        }
        else
        {
            const label patchi = mesh().boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(mesh().boundaryMesh()[patchi]))
            {
                label lpfI = fII - mesh().boundaryMesh()[patchi].start();
                const labelList& fcs = mesh().boundary()[patchi].faceCells();

                if (cRegion()[fcs[lpfI]] == 0)
                {
                    mPhiB[masterGIB_][fiBm] = mPhiB[patchi][lpfI];
                    fiBm += 1;
                }
                else
                {
                    mPhiB[slaveGIB_][fiBs] = mPhiB[patchi][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::makePhiTr() const
{
    if (phiTrPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "phiTr",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimVolume/dimTime
        );
    surfaceScalarField& phiTr = *phiTrPtr_;
    phiTr = mesh().phi() - phiB();
}


void Foam::fvMeshGIBChangersBase::makePhiTrS() const
{
    if (phiTrSPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrSPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "phiTrS",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh(),
            dimVolume/dimTime
        );
    surfaceScalarField& phiTrS = *phiTrSPtr_;
    scalarField& phiTrSInt = phiTrS.ref();
    surfaceScalarField::Boundary& phiTrSB = phiTrS.boundaryFieldRef();

    const boolList faceIndi0 = faceIndicator0();
    const boolList& popUpC = popUpCells();

    // Mark faces of the pop cells that they were not an interface
    // previously. These faces have to move and get shrinked volume
    boolList ppUp(mesh().points().size(), false);
    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            const labelList& cPoints = mesh().cellPoints()[celli];
            forAll(cPoints, cPointi)
            {
                ppUp[cPoints[cPointi]] = true;
            }
            const labelList& cFaces = mesh().cells()[celli];
            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                {
                    if (faceIndi0[facei] == true)
                    {
                        const labelList& f = mesh().faces()[facei];
                        forAll(f, pointi)
                        {
                            ppUp[f[pointi]] = false;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList(mesh(), ppUp, plusEqOp<bool>(), true);

    // Construction of the fluxes: based on the currect interface points and
    // the old interface points we make the pointfield with all the points
    // at the interface. By doint that we have shrinked volume.
    const boolList& interfacePoints = interP();
    const boolList& interfacePoints0 = interP0();

    vectorField oldPs = prevPoints_;
    vectorField newPs = prevPoints_;

    const vectorField& oldBl = recAllPoints0();

    boolList fPopIndi(mesh().faces().size(), false);
    forAll(mesh().points(), pointi)
    {
        const labelList& pFacesi = mesh().pointFaces()[pointi];
        if (ppUp[pointi])
        {
            forAll(pFacesi, pFacei)
            {
                const label facei = pFacesi[pFacei];
                if (facei < mesh().nInternalFaces())
                {
                    fPopIndi[facei] = true;
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchi])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchi])
                    )
                    {
                        fPopIndi[facei] = true;
                    }
                }
            }
        }
    }

    syncTools::syncFaceList(mesh(), fPopIndi, orEqOp<bool>());

    forAll(fPopIndi, facei)
    {
        if (fPopIndi[facei] == true)
        {
            const labelList& f = mesh().faces()[facei];
            forAll(f, fPointi)
            {
                const label pointi = f[fPointi];
                if
                (
                    (interfacePoints0[pointi] == false)
                 && (interfacePoints[pointi] == true)
                )
                {
                    newPs[pointi] = oldBl[pointi];
                }
            }
        }
    }

    tmp<scalarField> tsweptVols(new scalarField(mesh().faces().size(), 0));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(fPopIndi, facei)
    {
        if (fPopIndi[facei] == true)
        {
            sweptVols[facei] = mesh().faces()[facei].sweptVol(oldPs, newPs);
        }
    }

    dimensionedScalar deltaT = mesh().time().deltaT();
    phiTrSInt = scalarField::subField(sweptVols, mesh().nInternalFaces());
    phiTrSInt /= deltaT.value();

    const fvPatchList& patches = mesh().boundary();

    // Needs recheck
    forAll(patches, patchi)
    {
        if (!isA<indirectPolyPatch>(patches[patchi].patch()))
        {
            phiTrSB[patchi] = patches[patchi].patchSlice(sweptVols);
            phiTrSB[patchi] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, facei)
    {
        const label fII = fzNew[facei];
        if (fII < mesh().nInternalFaces())
        {
            phiTrSB[masterGIB_][facei] = phiTrSInt[fII];
            phiTrSB[slaveGIB_][facei] = -phiTrSInt[fII];
            if (fmNew[facei])
            {
                phiTrSB[masterGIB_][facei] *= -1;
                phiTrSB[slaveGIB_][facei] *= -1;
                fiBm += 1;
                fiBs += 1;
            }
        }
        else
        {
            const label patchi = mesh().boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(mesh().boundaryMesh()[patchi]))
            {
                label lpfI = fII - mesh().boundaryMesh()[patchi].start();
                const labelList& fcs = mesh().boundary()[patchi].faceCells();

                if (cRegion()[fcs[lpfI]] == 0)
                {
                    phiTrSB[masterGIB_][fiBm] = phiTrSB[patchi][lpfI];
                    fiBm += 1;
                }
                else
                {
                    phiTrSB[slaveGIB_][fiBs] = phiTrSB[patchi][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::makePhiTrG() const
{
    if (phiTrGPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    phiTrGPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "phiTrG",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh(),
            dimVolume/dimTime
        );
    surfaceScalarField& phiTrG = *phiTrGPtr_;
    scalarField& phiTrGInt = phiTrG.ref();
    surfaceScalarField::Boundary& phiTrGB = phiTrG.boundaryFieldRef();

    const boolList faceIndi = faceIndicator();

    const boolList& popUpC = popUpCells();

    // Mark faces of the pop cells that they were not an interface
    // previously. These faces have to move and get shrinked volume.
    boolList ppUp(mesh().points().size(), false);

    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            const labelList& cPoints = mesh().cellPoints()[celli];
            forAll(cPoints, cPointi)
            {
                ppUp[cPoints[cPointi]] = true;
            }
            const labelList& cFaces = mesh().cells()[celli];
            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                if (faceIndi[facei] == true)
                {
                    const labelList& f = mesh().faces()[facei];
                    forAll(f, fPointi)
                    {
                        ppUp[f[fPointi]] = false;
                    }
                }
            }
        }
    }

    syncTools::syncPointList(mesh(), ppUp, plusEqOp<bool>(), true);

    boolList fPopIndi(mesh().faces().size(), false);
    forAll(mesh().points(), pointi)
    {
        const labelList& pFacesi = mesh().pointFaces()[pointi];
        if (ppUp[pointi])
        {
            forAll(pFacesi, pFacei)
            {
                const label facei = pFacesi[pFacei];
                if (facei < mesh().nInternalFaces())
                {
                    fPopIndi[facei] = true;
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchi])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchi])
                    )
                    {
                        fPopIndi[facei] = true;
                    }
                }
            }
        }
    }

    syncTools::syncFaceList(mesh(), fPopIndi, orEqOp<bool>());

    // Construction of the fluxes: based on the currect interface points and
    // the old interface points we make the pointfield with all the points at
    // the interface. By doint that we have shrinked volume.
    //const boolList& interfacePoints = interP();
    const boolList& interfacePoints0 = interP0();

    vectorField oldPs = recAllPoints0();
    vectorField newPs = recAllPoints0();

    forAll(fPopIndi, facei)
    {
        if (fPopIndi[facei] == true)
        {
            const labelList& f = mesh().faces()[facei];
            forAll(f, fPointi)
            {
                const label pointi = f[fPointi];
                if (interfacePoints0[pointi] == true)
                {
                    newPs[pointi] = prevPoints_[pointi];
                }
            }
        }
    }

    tmp<scalarField> tsweptVols(new scalarField(mesh().faces().size(), 0));
    scalarField& sweptVols = tsweptVols.ref();

    forAll(fPopIndi, facei)
    {
        if (fPopIndi[facei] == true)
        {
            sweptVols[facei] = mesh().faces()[facei].sweptVol(newPs, oldPs);
        }
    }

    dimensionedScalar deltaT = mesh().time().deltaT();
    phiTrGInt = scalarField::subField(sweptVols, mesh().nInternalFaces());
    phiTrGInt /= deltaT.value();

    const fvPatchList& patches = mesh().boundary();

    // Needs recheck
    forAll(patches, patchi)
    {
        if (!isA<indirectPolyPatch>(patches[patchi].patch()))
        {
            phiTrGB[patchi] = patches[patchi].patchSlice(sweptVols);
            phiTrGB[patchi] /= deltaT.value();
        }
    }
    const labelList& fzNew = fl();
    const boolList& fmNew = fm();

    label fiBm = 0;
    label fiBs = 0;

    forAll(fmNew, facei)
    {
        const label fII = fzNew[facei];
        if (fII < mesh().nInternalFaces())
        {
            phiTrGB[masterGIB_][facei] = phiTrGInt[fII];
            phiTrGB[slaveGIB_][facei] = -phiTrGInt[fII];
            if (fmNew[facei])
            {
                phiTrGB[masterGIB_][facei] *= -1;
                phiTrGB[slaveGIB_][facei] *= -1;
            }
            fiBm += 1;
            fiBs += 1;
        }
        else
        {
            const label patchi = mesh().boundaryMesh().whichPatch(fII);
            if (!isA<indirectPolyPatch>(mesh().boundaryMesh()[patchi]))
            {
                const label lpfI = fII - mesh().boundaryMesh()[patchi].start();
                const labelList& fcs = mesh().boundary()[patchi].faceCells();

                if (cRegion()[fcs[lpfI]] == 0)
                {
                    phiTrGB[masterGIB_][fiBm] = phiTrGB[patchi][lpfI];
                    fiBm += 1;
                }
                else
                {
                    phiTrGB[slaveGIB_][fiBs] = phiTrGB[patchi][lpfI];
                    fiBs += 1;
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::makeFaceIndicator() const
{
    if (faceIndicatorPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    faceIndicatorPtr_ = new boolList(mesh().faces().size(), false);
    boolList& faceIndicator = *faceIndicatorPtr_;

    const labelList& fZone = fl();

    // Treatment of the pop up/in cells
    forAll(fZone, fZoneFacei)
    {
        const label facei = fZone[fZoneFacei];
        faceIndicator[facei] = true;
        if (!(facei < mesh().nInternalFaces()))
        {
            const label patchi = mesh().boundaryMesh().whichPatch(facei);
            if
            (
                isA<wedgeFvPatch>(mesh().boundary()[patchi])
             || isA<emptyFvPatch>(mesh().boundary()[patchi])
            )
            {
                faceIndicator[facei] = false;
            }
        }
    }

    syncTools::syncFaceList(mesh(), faceIndicator, orEqOp<bool>());
}


void Foam::fvMeshGIBChangersBase::makeFaceIndicator0() const
{
    if (faceIndicator0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    faceIndicator0Ptr_ = new boolList(mesh().faces().size(), false);
    boolList& faceIndicator0 = *faceIndicator0Ptr_;

    const labelList& fZone0 = fl0();
    forAll(fZone0, fZone0Facei)
    {
        const label facei = fZone0[fZone0Facei];
        faceIndicator0[facei] = true;
        if (!(facei < mesh().nInternalFaces()))
        {
            const label patchi = mesh().boundaryMesh().whichPatch(facei);
            if (!mesh().boundary()[patchi].coupled())
            {
                faceIndicator0[facei] = false;
            }
        }
    }

    syncTools::syncFaceList(mesh(), faceIndicator0, orEqOp<bool>());
}


void Foam::fvMeshGIBChangersBase::makecRegion() const
{
    if (cRegionPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList blockedFace(mesh().nFaces(), false);

    const labelList& fList = fl();
    forAll(fList, facei)
    {
        blockedFace[fList[facei]] = true;
    }

    // Parallel sync !imporant for regionSplit
    syncTools::syncFaceList(mesh(), blockedFace, orEqOp<bool>());

    regionSplit cellMark(mesh(), blockedFace);
    labelList cellIndi(cellMark);
    modifyRegionLabels(cellIndi);
    cRegionPtr_ = new labelList(cellIndi);
}


void Foam::fvMeshGIBChangersBase::makecRegion0() const
{
    if (cRegion0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boolList blockedFace(mesh().nFaces(), false);

    const labelList& fList = fl0();
    forAll(fList, facei)
    {
        blockedFace[fList[facei]] = true;
    }

    // Parallel sync !imporant for regionSplit
    syncTools::syncFaceList(mesh(), blockedFace, orEqOp<bool>());

    regionSplit cellMark(mesh(), blockedFace);
    labelList cellIndi(cellMark);
    modifyRegionLabels(cellIndi);
    cRegion0Ptr_ = new labelList(cellIndi);
}


void Foam::fvMeshGIBChangersBase::makeFullSnapCells() const
{
    if (fullSnapCellsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    const boolList& interfacePoints = interP();

    DynamicList<label> dlc(mesh().cells().size());

    forAll(mesh().cells(), celli)
    {
        const labelList& cPoints = mesh().cellPoints()[celli];
        bool foundUnsnapped = false;
        forAll(cPoints, cPointi)
        {
            if (!interfacePoints[cPoints[cPointi]])
            {
                foundUnsnapped = true;
            }
        }
        if (!foundUnsnapped)
        {
            dlc.append(celli);
        }
    }
    dlc.shrink();

    fullSnapCellsPtr_ = new labelList(dlc);
}


void Foam::fvMeshGIBChangersBase::makeFullSnapCellPoints() const
{
    if (fullSnapCellPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    const labelList& fsc = fullSnapCells();
    boolList markPoints(mesh().points().size(), false);
    forAll(fsc, celli)
    {
        const labelList& cPoints = mesh().cellPoints()[fsc[celli]];
        forAll(cPoints, cPointi)
        {
            markPoints[cPoints[cPointi]] = true;
        }
    }

    syncTools::syncPointList(mesh(), markPoints, plusEqOp<bool>(), true);

    DynamicList<label> dlp(mesh().points().size());
    forAll(markPoints, pointi)
    {
        if (markPoints[pointi])
        {
            dlp.append(pointi);
        }
    }

    dlp.shrink();

    fullSnapCellPointsPtr_ = new labelList(dlp);
}


void Foam::fvMeshGIBChangersBase::makePopCellPoints() const
{
    if (popCellPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popCellPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& popCellPoints = *popCellPointsPtr_;

    const boolList& popUpC = popUpCells();

    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            const labelList& cPoints = mesh().cellPoints()[celli];
            forAll(cPoints, cPointi)
            {
                popCellPoints[cPoints[cPointi]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh(), popCellPoints, plusEqOp<bool>(), true);
}


void Foam::fvMeshGIBChangersBase::makePopSPoints() const
{
    if (popSPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popSPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& popSPoints = *popSPointsPtr_;

    const boolList faceIndi0 = faceIndicator0();

    const boolList& popUpC = popUpCells();

    // Mark points at the pop cells and on the interface in the previous
    // iteration
    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            const labelList& cPoints = mesh().cellPoints()[celli];
            forAll(cPoints, cPointi)
            {
                popSPoints[cPoints[cPointi]] = true;
            }
            const labelList& cFaces = mesh().cells()[celli];
            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                if (faceIndi0[facei] == true)
                {
                    const labelList& f = mesh().faces()[facei];
                    forAll(f, fPointi)
                    {
                        popSPoints[f[fPointi]] = false;
                    }
                }
            }
        }
    }

    syncTools::syncPointList(mesh(), popSPoints, plusEqOp<bool>(), true);
}


void Foam::fvMeshGIBChangersBase::makePopGPoints() const
{
    if (popGPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    popGPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& popGPoints = *popGPointsPtr_;

    const boolList faceIndi = faceIndicator();

    const boolList& popUpC = popUpCells();

    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            const labelList& cPoints = mesh().cellPoints()[celli];
            forAll(cPoints, cPointi)
            {
                popGPoints[cPoints[cPointi]] = true;
            }
            const labelList& cFaces = mesh().cells()[celli];
            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                if (faceIndi[facei] == true)
                {
                    const labelList& f = mesh().faces()[facei];
                    forAll(f, fPointi)
                    {
                        popGPoints[f[fPointi]] = false;
                    }
                }
            }
        }
    }

    syncTools::syncPointList(mesh(), popGPoints, plusEqOp<bool>(), true);
}


void Foam::fvMeshGIBChangersBase::makeRecPoints0() const
{
    if (recPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recPoints0 = *recPoints0Ptr_;

    tmp<pointField> tsurP(triSM().points());
    const pointField& surP = tsurP();
    const pointField& surP0 = triSM().triSurface::points0();
    const labelList& hitInd = hitIndex();
    const labelList& addrM = mesh().boundary()[masterId()].patch().meshPoints();

    boolList oldnewintPoints(mesh().points().size(), false);
    boolList marknewFaces(mesh().faces().size(), false);
    boolList markoldFaces(mesh().faces().size(), false);

    const labelList& fZone0 = fl0();
    const labelList& fZone = fl();
    forAll(fZone0, facei)
    {
        markoldFaces[fZone0[facei]] = true;
    }
    forAll(fZone, facei)
    {
        marknewFaces[fZone[facei]] = true;
    }
    forAll(mesh().faces(), facei)
    {
        if (marknewFaces[facei] && markoldFaces[facei])
        {
            forAll(mesh().faces()[facei], fPointi)
            {
                const label pointi = mesh().faces()[facei][fPointi];
                oldnewintPoints[pointi] = true;
            }
        }
    }

    // We want to express the point based on the coordinates of the
    // points of the hit triangle of the stl.
    // p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a 3x3 system
    forAll(addrM, pI)
    {
        label addrI = addrM[pI];
        if (!oldnewintPoints[addrI])
        {
            label hitIndexI = hitInd[addrI];
            const triFace& triList =
                triSM().triSurface::operator[](hitIndexI);

            vector p = hitPoint()[addrI];
            vector p0_1 = surP0[triList[0]];
            vector p0_2 = surP0[triList[1]];
            vector p0_3 = surP0[triList[2]];

            vector w = triList.findTriangleWeights(p, surP);

            recPoints0[addrI] = w[0]*p0_1 + w[1]*p0_2 + w[2]*p0_3;
        }
        else
        {
            recPoints0[addrI] = prevPoints_[addrI];
        }
    }

    syncPoints(recPoints0);
}


void Foam::fvMeshGIBChangersBase::makeRecAllPoints0() const
{
    if (recAllPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recAllPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recAllPoints0 = *recAllPoints0Ptr_;

    pointField disp (mesh().points().size(), Zero);

    if (true)
    {
        tmp<pointField> tsurP(triSM().points());
        const pointField& surP = tsurP();
        const pointField& surP0 = triSM().triSurface::points0();
        const labelList& hitInd = hitIndex();
        const labelList& addr =
            mesh().boundary()[masterId()].patch().meshPoints();

        if (false)
        {
            const word timeName(mesh().time().timeName());
            simpleVTKWriter(surP).write("stl_" + timeName + ".vtk");
            simpleVTKWriter(surP0).write("stl0_" + timeName + ".vtk");
        }

        // We want to express the point based on the coordinates of the
        // points of the hit triangle of the stl.
        // p = w1*p1 + w2*p2 + w3*p3  -->  we are solving a 3x3 system.
        // First for the master points.
        boolList pointMoved(mesh().points().size(), false);

        forAll(addr, pI)
        {
            label addrI = addr[pI];
            label hitIndexI = hitInd[addrI];
            const triFace& triList = triSM().triSurface::operator[](hitIndexI);
            vector p = hitPoint()[addrI];
            vector p0_0 = surP0[triList[0]];
            vector p0_1 = surP0[triList[1]];
            vector p0_2 = surP0[triList[2]];

            vector w = triList.findTriangleWeights(p, surP);
            vector newLoc = w[0]*p0_0 + w[1]*p0_1 + w[2]*p0_2;
            disp[addrI] = newLoc - recAllPoints0[addrI];
            pointMoved[addrI] = true;
        }

        syncTools::syncPointList(mesh(), pointMoved, plusEqOp<bool>(), true);

        const labelList& addrS =
            mesh().boundary()[slaveId()].patch().meshPoints();

        // Second the points in slave patch that they have not been updated
        // from the master loop points
        forAll(addrS, pI)
        {
            label addrI = addrS[pI];
            if (pointMoved[addrI] == false)
            {
                label hitIndexI = hitInd[addrI];
                const triFace& triList = triSM().triSurface::operator[](hitIndexI);
                vector p = hitPoint()[addrI];
                vector p0_1 = surP0[triList[0]];
                vector p0_2 = surP0[triList[1]];
                vector p0_3 = surP0[triList[2]];

                vector w = triList.findTriangleWeights(p, surP);

                vector newLoc = w[0]*p0_1 + w[1]*p0_2 + w[2]*p0_3;
                disp[addrI] = newLoc - recAllPoints0[addrI];
                pointMoved[addrI] = true;
            }
        }

        const labelList& fscp = fullSnapCellPoints();

        syncTools::syncPointList(mesh(), disp, maxMagEqOp(), vector::zero);
        recAllPoints0 += disp;

        const boolList& interPoints = interP();
        const boolList& interPoints0 = interP0();

        computeOldPositionsInUnsnappedCells(recAllPoints0, fscp);

        const boolList& meshQualityP = meshQualityPoints();

        forAll(meshQualityP, pointi)
        {
            if
            (
                meshQualityP[pointi]
             && interPoints[pointi]
             && interPoints0[pointi]
            )
            {
                recAllPoints0[pointi] = prevPoints_[pointi];
            }
        }
    }
    else
    {
        recAllPoints0 = prevPoints_;
    }
}


void Foam::fvMeshGIBChangersBase::makePopUpCells() const
{
    const labelList& cReg = cRegion();
    const labelList& cReg0 = cRegion0();

    popUpCellsPtr_ = new boolList(mesh().cells().size(), false);
    boolList& popUpCells = *popUpCellsPtr_;

    forAll(mesh().cells(), celli)
    {
        if (cReg[celli] != cReg0[celli])
        {
            popUpCells[celli] = 1;
        }
        else
        {
            popUpCells[celli] = 0;
        }
    }
    if (mesh().time().outputTime() && debugMode_)
    {
        writeScalarField("popCells", popUpCells);
    }
}


void Foam::fvMeshGIBChangersBase::makeBoundaryPoints() const
{
    if (boundaryPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    boundaryPointsPtr_ = new boolList(mesh().points().size(), false);

    boolList& boundaryPoints = *boundaryPointsPtr_;

    for (label facei = mesh().nInternalFaces(); facei < mesh().nFaces(); facei++)
    {
        const label patchi = mesh().boundaryMesh().whichPatch(facei);
        if
        (
            !isA<indirectPolyPatch>(mesh().boundaryMesh()[patchi])
         && !isA<wedgeFvPatch>(mesh().boundary()[patchi])
         && !isA<emptyFvPatch>(mesh().boundary()[patchi])
         && !mesh().boundary()[patchi].coupled()
        )
        {
            const labelList& f = mesh().faces()[facei];
            forAll(f, fPointi)
            {
                boundaryPoints[f[fPointi]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh(), boundaryPoints, plusEqOp<bool>(), true);
}


void Foam::fvMeshGIBChangersBase::extendBoundaryBaseMesh
(
    const labelList& interestFaces,
    vectorField& neiCC
) const
{
    const vectorField& baseCf = *baseCf_;
    const vectorField& baseCC = *baseCC_;
    forAll(interestFaces, facei)
    {
        label interestFacesI = interestFaces[facei];
        if (interestFacesI >= mesh().nInternalFaces())
        {
            const label patchi =
                mesh().boundaryMesh().whichPatch(interestFacesI);
            const polyPatch& pp = mesh().boundary()[patchi].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                if (!pp.coupled())
                {
                    const label own = mesh().faceOwner()[interestFacesI];
                    const vector& ccI = baseCC[own];

                    vector vfc = baseCf[interestFacesI] - ccI;
                    scalar mvfc = mag(vfc);
                    const vector& nf = (vfc)/(mvfc+SMALL);

                    // Extend base mesh face center points by cfCC vector
                    neiCC[facei] = baseCf[interestFacesI]+(nf*mvfc);
                }
            }
        }
    }
    if (false)
    {
        const word timeName(mesh().time().timeName());
        simpleVTKWriter(neiCC).write("extBoundary_" + timeName + ".vtk");
    }
}


Foam::List<Foam::face> Foam::fvMeshGIBChangersBase::preparePatch
(
    const fvPatch& gibPatch
)
{
    faceList faces(gibPatch.patch().localFaces());

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    const boolList& flipmap = gibPolyPatch.fm();

    forAll(faces, facei)
    {
        if (flipmap[facei])
        {
            faces[facei].flip();
        }
    }

    return faces;
}


void Foam::fvMeshGIBChangersBase::syncPoints(pointField& points) const
{
    // Syncing the processor boundary points
    pointField np1(points);
    pointField np2(points);

    syncTools::syncPointPositions
    (
        mesh(),
        np1,
        minEqOp<point>(),             // combine op
        point(GREAT, GREAT, GREAT)    // null
    );
    syncTools::syncPointPositions
    (
        mesh(),
        np2,
        maxEqOp<point>(),                // combine op
        point(-GREAT, -GREAT, -GREAT)    // null
    );

    forAll(points, pointi)
    {
        points[pointi] = (np1[pointi] + np2[pointi])/2;
    }
}


Foam::List<Foam::label>
Foam::fvMeshGIBChangersBase::checkingIntersectingFaces() const
{
    const labelList& fl = fl0();

    boolList interestedFace(mesh().nFaces(), false);
    boolList constraintFace(mesh().nFaces(), false);

    // If initiallization take all the faces
    label ggibSize = fl.size();

    reduce(ggibSize, sumOp<label>());

    // Mark all the faces we dont want to be included
    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& pp = mesh().boundary()[patchi].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if
            (
                isA<wedgeFvPatch>(mesh().boundary()[patchi])
             || isA<emptyFvPatch>(mesh().boundary()[patchi])
            )
            {
                forAll(pp, pFacei)
                {
                    const label facei = pp.start() + pFacei;
                    constraintFace[facei] = true;
                }
            }
        }
    }

    if (ggibSize == 0)
    {
        // If old time doens't exist (for example when runing createGIB)
        // add all the faces except the constraint faces
        forAll(interestedFace, facei)
        {
            if (constraintFace[facei] == false)
            {
                interestedFace[facei] = true;
            }
        }
        forAll(mesh().boundary(), patchi)
        {
            const polyPatch& pp = mesh().boundary()[patchi].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                if ((!pp.coupled()) && !includeWalls())
                {
                    forAll(pp, pFacei)
                    {
                        const label facei = pp.start() + pFacei;
                        interestedFace[facei] = false;
                    }
                }
            }
        }
        DynamicList<label> dfl(mesh().faces().size());

        forAll(interestedFace, facei)
        {
            if (interestedFace[facei])
            {
                dfl.append(facei);
            }
        }
        dfl.shrink();
        return labelList(dfl);
    }

    // Mark all the faces of the previous interface
    forAll(fl, facei)
    {
        interestedFace[fl[facei]] = true;
    }

    // Mark all the boundaries
    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& pp = mesh().boundary()[patchi].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            forAll(pp, pFacei)
            {
                const label facei = pp.start() + pFacei;
                if (!constraintFace[facei])
                {
                    interestedFace[facei] = true;
                }
            }
        }
    }

    boolList bandFaces(mesh().nFaces(), false);

    // Controling the band
    for (int i = 1; i < 4; i++)
    {
        // Mark all the face-cell-faces of the internal faces
        forAll(mesh().neighbour(), facei)
        {
            if (interestedFace[facei])
            {
                const label own = mesh().owner()[facei];
                const label nei = mesh().neighbour()[facei];
                {
                    const labelList& cFaces = mesh().cells()[own];
                    forAll(cFaces, cFacei)
                    {
                        const label cfII = cFaces[cFacei];
                        if (!constraintFace[cfII])
                        {
                            bandFaces[cfII] = true;
                        }
                    }
                }
                {
                    const labelList& cFaces = mesh().cells()[nei];
                    forAll(cFaces, cFacei)
                    {
                        const label cfII = cFaces[cFacei];
                        if (!constraintFace[cfII])
                        {
                            bandFaces[cfII] = true;
                        }
                    }
                }
            }
        }

        // Mark all the face-cell-faces of the boundaryCells
        forAll(mesh().boundary(), patchi)
        {
            const polyPatch& pp = mesh().boundary()[patchi].patch();
            if (!isA<indirectPolyPatch>(pp))
            {
                const labelList& fcs = mesh().boundary()[patchi].faceCells();
                forAll(pp, pFacei)
                {
                    const label facei = pp.start() + pFacei;
                    if (interestedFace[facei])
                    {
                        const label fCelli = fcs[pFacei];
                        const labelList& cFaces = mesh().cells()[fCelli];
                        forAll(cFaces, cFacei)
                        {
                            const label cfII = cFaces[cFacei];
                            if (!constraintFace[cfII])
                            {
                                bandFaces[cfII] = true;
                            }
                        }
                    }
                }
            }
        }

        forAll(interestedFace, facei)
        {
            if (bandFaces[facei])
            {
                interestedFace[facei] = true;
            }
        }
    }

    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& pp = mesh().boundary()[patchi].patch();
        if (!isA<indirectPolyPatch>(pp))
        {
            if ((!pp.coupled()) && !includeWalls())
            {
                forAll(pp, pFacei)
                {
                    const label facei = pp.start() + pFacei;
                    interestedFace[facei] = false;
                }
            }
        }
    }

    DynamicList<label> dfl(mesh().faces().size());

    forAll(interestedFace, facei)
    {
        if (interestedFace[facei])
        {
            dfl.append(facei);
        }
    }
    dfl.shrink();

    return labelList(dfl);
}


bool Foam::fvMeshGIBChangersBase::checkIfEdgeToleranceIntersecting
(
    const face& face1,
    const face& face2
) const
{
    int samePoints = 0;
    bool tolIssue = false;

    forAll(face1, f1)
    {
        forAll(face2, f2)
        {
            if (face1[f1] == face2[f2])
            {
                samePoints +=1;
            }
        }
    }
    if (samePoints == 2)
    {
        tolIssue = true;
    }

    return tolIssue;
}


void Foam::fvMeshGIBChangersBase::doTangentialBoundaryMotion
(
    pointField& newPoints
) const
{
    const vectorField& baseSf = *baseSf_;
    const pointField& basePoints = *basePoints_;

    const boolList& bPoints = boundaryPoints();

    const labelList& fList = fl();

    boolList gibPoints(basePoints.size(), false);

    forAll(fList, facei)
    {
        const label fListi = fList[facei];
        const face& f = mesh().faces()[fListi];
        forAll(f, fPointi)
        {
            gibPoints[f[fPointi]] = true;
        }
    }

    boolList gibBoundaryPoints(basePoints.size(), false);

    forAll(gibBoundaryPoints, pointi)
    {
        if (bPoints[pointi] && gibPoints[pointi])
        {
            gibBoundaryPoints[pointi] = true;
        }
    }

    const labelListList& pFaces = mesh().pointFaces();
    forAll(gibBoundaryPoints, pointi)
    {
        if (gibBoundaryPoints[pointi])
        {
            vector dis = newPoints[pointi] - basePoints[pointi];

            const labelList& pFacesi = pFaces[pointi];
            forAll(pFacesi, pFacei)
            {
                const label facei = pFacesi[pFacei];
                if (facei >= mesh().nInternalFaces())
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    const polyPatch& poly = mesh().boundaryMesh()[patchi];
                    if (isA<directPolyPatch>(poly))
                    {
                        if
                        (
                           !(
                                mesh().boundaryMesh()[patchi].coupled()
                             || isA<wedgePolyPatch>(mesh().boundaryMesh()[patchi])
                             || isA<emptyPolyPatch>(mesh().boundaryMesh()[patchi])
                            )
                        )
                        {
                            vector n = baseSf[facei]/mag(baseSf[facei]);
                            dis -= (dis&n)*n;
                        }
                    }
                }
            }
            newPoints[pointi] = basePoints[pointi] + dis;
        }
    }
}


Foam::tmp<Foam::pointField> Foam::fvMeshGIBChangersBase::findSnappedPoints
(
    const bool& fromBase
) const
{
    const pointField& baseP(*basePoints_);

    // Initialize the snapped with the current points
    tmp<pointField> npt(new pointField(mesh().points()));
    pointField& np = npt.ref();
    if (fromBase)
    {
        np = baseP;
    }
    else
    {
        np = mesh().points();
    }

    const boolList& interPoints = interP();
    const boolList& interPoints0 = interP0();
    DynamicList<label> dIntPoints(mesh().points().size());

    forAll(interPoints, pointi)
    {
        if (interPoints[pointi])
        {
            dIntPoints.append(pointi);
        }
    }
    dIntPoints.shrink();

    if (!fromBase)
    {
        // The points of the current GIB are reverted to the basePoint location
        forAll(dIntPoints, pointi)
        {
            const label gpI = dIntPoints[pointi];
            np[gpI] = baseP[gpI];
        }
    }
    pointField dIntPointLoc = pointField(dIntPoints.size(), Zero);

    List<pointIndexHit> nearest;

    boolList multInterPoints = findMultInterPoints();

    const boolList& popCellP = popCellPoints();

    // Mark points which belong to the cells that all their points are snapped
    const labelList& fscp = fullSnapCellPoints();
    boolList gfullSnapCellPoints(mesh().points().size(), false);
    forAll(fscp, pointi)
    {
        gfullSnapCellPoints[fscp[pointi]] = true;
    }

    const boolList& bPoints = boundaryPoints();

    forAll(dIntPoints, pointi)
    {
        dIntPointLoc[pointi] = baseP[dIntPoints[pointi]];
    }

    smoothInternalBasePoints(dIntPointLoc, dIntPoints);

    forAll(dIntPoints, pointi)
    {
        label gpI = dIntPoints[pointi];
        if
        (
            (
                (
                    !interPoints0[gpI]
                 || multInterPoints[gpI]
                 || gfullSnapCellPoints[gpI]
                 || popCellP[gpI]
                )
            && !bPoints[gpI]
            )
         || fromBase
        )
        {}
        else
        {
            dIntPointLoc[pointi] = mesh().points()[gpI];
        }
    }

    triSM().findNearest
    (
        dIntPointLoc,
        scalarField(dIntPointLoc.size(), sqr(GREAT)),
        nearest
    );

    labelIOList& hitIndex = *hitIndexPtr_;
    hitIndex.resize(mesh().points().size(), -1);

    hitPointPtr_ = new vectorField(mesh().points().size(), Zero);
    vectorField& hitPoint = *hitPointPtr_;

    forAll(dIntPoints, pointi)
    {
        const pointIndexHit& nearI = nearest[pointi];
        np[dIntPoints[pointi]] = nearI.hitPoint();
        hitIndex[dIntPoints[pointi]] = nearI.index();
    }

    syncPoints(np);
    hitPoint = np;

    pointField fromPoints(mesh().points().size(), Zero);
    forAll(dIntPoints, pointi)
    {
        const label gpI = dIntPoints[pointi];
        fromPoints[gpI] = dIntPointLoc[pointi];
    }

    //checkConcaveBoundaryPoints(baseP, np);

    fullSnappedPointsTreatment(np, baseP);

    unsnapForQuality(np, baseP);

    syncPoints(np);

    return npt;
}


void Foam::fvMeshGIBChangersBase::unsnapForQuality
(
    pointField& newPoints,
    const pointField& fromPoints
) const
{
    pointField snappedPoints = newPoints;
    const boolList& isInterfacePoint = interP();
    scalarField r(newPoints.size(), 0);
    const label maxUnsnapIntervals = 10;
    label unsnapStep;

    if (meshQualityPointsPtr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }
    meshQualityPointsPtr_ = new boolList(mesh().points().size(), false);
    boolList& meshQualityPoints = *meshQualityPointsPtr_;

    boolList hasInterPoint(mesh().faces().size(), false);
    forAll(hasInterPoint, facei)
    {
        const face& f = mesh().faces()[facei];
        forAll(f, fPointi)
        {
            if (isInterfacePoint[f[fPointi]])
            {
                hasInterPoint[facei] = true;
                break;
            }
        }
    }

    for (unsnapStep = 0; unsnapStep < maxUnsnapIntervals; unsnapStep++)
    {
        bool qualityOK = true;
        scalarList nonOrth(mesh().faces().size(), GREAT);
        scalarList ownPyrVol(mesh().faces().size(), scalar(0));
        scalarList neiPyrVol(mesh().faces().size(), scalar(0));

        vectorField d(mesh().faces().size(), Zero);
        forAll(mesh().faces(), facei)
        {
            if (hasInterPoint[facei])
            {
                label celli = mesh().faceOwner()[facei];
                const face& f = mesh().faces()[facei];
                point ownCentre
                (
                    mesh().cells()[celli].centre(newPoints, mesh().faces())
                );

                // Calc cell-centre-to-cell-centre vectors
                d[facei] = -ownCentre;

                // Calculate pyramid volumes
                point fc(f.centre(newPoints));
                vector fa(f.areaNormal(newPoints));
                ownPyrVol[facei] =
                   -primitiveMeshTools::pyramidVol(ownCentre, fc, fa);

                if (facei < mesh().faceNeighbour().size())
                {
                    label cellj = mesh().faceNeighbour()[facei];
                    point neiCentre
                    (
                        mesh().cells()[cellj].centre(newPoints, mesh().faces())
                    );

                    d[facei] += neiCentre;

                    neiPyrVol[facei] =
                        primitiveMeshTools::pyramidVol(neiCentre, fc, fa);
                }
                else
                {
                    // This gets cancelled for coupled boundaries
                    d[facei] += f.centre(newPoints);
                }
            }
        }

        // Calculate across coupled boundaries
        syncTools::syncFaceList(mesh(), d, minusEqOp<vector>());

        forAll(mesh().faces(), facei)
        {
            if (hasInterPoint[facei])
            {
                // Calc non-orth
                const face& f = mesh().faces()[facei];
                vector nf =
                    f.areaNormal(newPoints)/stabilise(f.mag(newPoints), SMALL);
                nonOrth[facei] = (nf&d[facei])/stabilise(mag(d[facei]), SMALL);
            }
        }

        boolList unsnapP(newPoints.size(), false);
        forAll(mesh().faces(), facei)
        {
            if (hasInterPoint[facei])
            {
                if
                (
                    nonOrth[facei] < maxNonOrth_
                 || ownPyrVol[facei] < 0
                 || neiPyrVol[facei] < 0
                )
                {
                    qualityOK = false;
                    const face& f = mesh().faces()[facei];
                    forAll(f, fPointi)
                    {
                        unsnapP[f[fPointi]] = true;
                    }
                }
            }
        }

        reduce(qualityOK, andOp<bool>());
        if (qualityOK)
        {
            break;
        }

        // Sync here as unsnapP might have got changed by a non-coupled face
        // on one side
        syncTools::syncPointList(mesh(), unsnapP, orEqOp<bool>(), false);
        forAll(meshQualityPoints, pointi)
        {
            meshQualityPoints[pointi] =
                meshQualityPoints[pointi] || unsnapP[pointi];
        }

        forAll(unsnapP, pointi)
        {
            if (unsnapP[pointi])
            {
                r[pointi] = scalar(unsnapStep+1)/(maxUnsnapIntervals);
                newPoints[pointi] =
                    (
                        snappedPoints[pointi]
                      - (snappedPoints[pointi] - fromPoints[pointi])*r[pointi]
                    );
            }
        }
    }
    if (debugMode_)
    {
        scalar totalr(0), maxr(0);
        label numIntP(0);
        forAll(r, pointi)
        {
            if (isInterfacePoint[pointi])
            {
                totalr += r[pointi];
                maxr = max(r[pointi], maxr);
                numIntP++;
            }
        }
        Foam::reduce
        (
            std::tie(totalr, maxr, numIntP),
            ParallelOp<sumOp<scalar>, maxOp<scalar>, sumOp<label>>{},
            mesh().comm()
        );

        if (!numIntP)
        {
            // Avoid /0 below
            numIntP = 1;
        }
        Info<< "GIB mesh-quality constraint: Unsnapping fraction  "
            << "Mean = " << totalr/numIntP << ", max = " << maxr << endl;
    }
}


void Foam::fvMeshGIBChangersBase::fullSnappedPointsTreatment
(
    pointField& newPoints,
    const pointField& fromPoints
) const
{
    const labelList& fscp = fullSnapCellPoints();

    forAll(fscp, pointi)
    {
        const label gpI = fscp[pointi];
        newPoints[gpI] -= unsnapVar_*(newPoints[gpI] - fromPoints[gpI]);
    }
}


void Foam::fvMeshGIBChangersBase::computeOldPositionsInUnsnappedCells
(
    pointField& recAllPoints0,
    const labelList& fscp
) const
{
    forAll(fscp, pointi)
    {
        const label gpI = fscp[pointi];
        recAllPoints0[gpI] = mesh().points()[gpI];
    }
}


Foam::label Foam::fvMeshGIBChangersBase::gibFaceZone() const
{
    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>
        (
            mesh().boundary()[masterId()].patch()
        );

    return gibPolyPatch.zoneId();
}


Foam::label Foam::fvMeshGIBChangersBase::findOrCreateCellZone() const
{
    const label fZoneId = gibFaceZone();
    const word faceZoneName = mesh().faceZones()[fZoneId].name();

    const word cellZoneName = "inactive_"+faceZoneName;
    const label cZoneId = mesh().cellZones().findZoneID(cellZoneName);

    if (cZoneId != -1)
    {
        return cZoneId;
    }
    else
    {
        Info<< "cellZone was not found:" << endl;
        Info<< "Constructing new cellZone called " << cellZoneName << endl;
        const polyMesh& pMesh = mesh();
        polyMesh& cpMesh = const_cast<polyMesh&>(pMesh);
        label oldZoneSize = cpMesh.cellZones().size();
        cpMesh.cellZones().setSize(oldZoneSize + 1);

        label index = oldZoneSize;
        cpMesh.cellZones().set
        (
            index,
            new cellZone
            (
                cellZoneName,
                labelUList(),
                index,
                cpMesh.cellZones()
            )
        );

        return index;
    }
}


void Foam::fvMeshGIBChangersBase::calculateFlipMap() const
{
    updateSolidCellZone();
    const cellZone& cZone = mesh().cellZones()[findOrCreateCellZone()];
    boolList markCells(mesh().cells().size(), false);
    const labelList& fList = fl();
    boolList fm(fList.size(), false);

    forAll(cZone, celli)
    {
        markCells[cZone[celli]] = true;
    }

    forAll(fList, facei)
    {
        const label own = mesh().faceOwner()[fList[facei]];
        if (markCells[own])
        {
            fm[facei] = true;
        }
    }
    fmPtr_ = new boolList(fm);
}


void Foam::fvMeshGIBChangersBase::regionVisDebug() const
{
    const labelList& cReg = cRegion();

    volScalarField cellRegionF
    (
        IOobject
        (
            "cellRegionF",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimless, 0),
        "calculated"
    );

    forAll(cellRegionF, celli)
    {
        cellRegionF[celli] = cReg[celli];
    }

    if (mesh().time().outputTime())
    {
        cellRegionF.write();
    }
}


Foam::List<bool> Foam::fvMeshGIBChangersBase::findMultInterPoints() const
{
    boolList multInterPoints(mesh().points().size(), false);
    const boolList& multInterFaces = *multInterFacesPtr_;
    forAll(multInterFaces, facei)
    {
        if (multInterFaces[facei] == true)
        {
            const labelList& f = mesh().faces()[facei];
            forAll(f, fPointi)
            {
                multInterPoints[f[fPointi]] = true;
            }
        }
    }

    return multInterPoints;
}


void Foam::fvMeshGIBChangersBase::popShrinkFields()
{
    const boolList& popSP = popSPoints();

    // Mark all the faces and cells that are going to be popped
    // from the poped points
    boolList fPopIndi(mesh().faces().size(), false);
    boolList popIndi(mesh().cells().size(), false);
    forAll(mesh().points(), pointi)
    {
        const labelList& pFacesi = mesh().pointFaces()[pointi];
        const labelList& pCellsI = mesh().pointCells()[pointi];
        if (popSP[pointi])
        {
            forAll(pFacesi, pFacei)
            {
                const label facei = pFacesi[pFacei];
                if (facei < mesh().nInternalFaces())
                {
                    fPopIndi[facei] = true;
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchi])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchi])
                    )
                    {
                        fPopIndi[facei] = true;
                    }
                }
            }
            forAll(pCellsI, celli)
            {
                popIndi[pCellsI[celli]] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh(), fPopIndi, orEqOp<bool>());

    // Add actual face lables in the fPop list from marked boolList
    DynamicList<label> fPopD(mesh().faces().size());
    forAll(mesh().faces(), facei)
    {
        if (fPopIndi[facei])
        {
            fPopD.append(facei);
        }
    }
    fPopD.shrink();
    labelList fPop(fPopD);

    const surfaceScalarField& popShrinkPhi = phiTrS();

    shrinkPopFields<scalar, fvPatchField, volMesh>
    (
        fPop,
        popShrinkPhi,
        popIndi
    );
    shrinkPopFields<vector, fvPatchField, volMesh>
    (
        fPop,
        popShrinkPhi,
        popIndi
    );
    shrinkPopFields<sphericalTensor, fvPatchField, volMesh>
    (
        fPop,
        popShrinkPhi,
        popIndi
    );
    shrinkPopFields<symmTensor, fvPatchField, volMesh>
    (
        fPop,
        popShrinkPhi,
        popIndi
    );
    shrinkPopFields<tensor, fvPatchField, volMesh>
    (
        fPop,
        popShrinkPhi,
        popIndi
    );
}


void Foam::fvMeshGIBChangersBase::popGrowFields()
{
    const scalarField& cV0 = mesh().V0();

    const boolList& popGP = popGPoints();

    const boolList& popUpC = popUpCells();

    boolList fPopIndi(mesh().faces().size(), false);
    boolList popIndi(mesh().cells().size(), false);
    forAll(mesh().points(), pointi)
    {
        const labelList& pFacesi = mesh().pointFaces()[pointi];
        const labelList& pCellsI = mesh().pointCells()[pointi];
        if (popGP[pointi])
        {
            forAll(pFacesi, pFacei)
            {
                const label facei = pFacesi[pFacei];
                if (facei < mesh().nInternalFaces())
                {
                    fPopIndi[facei] = true;
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchi])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchi])
                    )
                    {
                        fPopIndi[facei] = true;
                    }
                }
            }
            forAll(pCellsI, celli)
            {
                popIndi[pCellsI[celli]] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh(), fPopIndi, orEqOp<bool>());

    DynamicList<label> fPopD(mesh().faces().size());
    forAll(mesh().faces(), facei)
    {
        if (fPopIndi[facei])
        {
            fPopD.append(facei);
        }
    }
    fPopD.shrink();
    labelList fPop(fPopD);

    const surfaceScalarField& phiPop = phiTrG();

    scalarField sV(mesh().cells().size(), 0.0);

    // oldVgrow -> volume when we start growing from 0 sized volume
    scalarField oldVolgrow = cV0;

    {
        forAll(fPop, facei)
        {
            const label fPopI = fPop[facei];
            if (fPopI < mesh().nInternalFaces())
            {
                scalar sVf = phiPop[fPopI]*mesh().time().deltaTValue();

                const label own = mesh().owner()[fPopI];
                const label nei = mesh().neighbour()[fPopI];
                sV[own] += sVf;
                sV[nei] -= sVf;
            }
            else
            {
                const label own = mesh().faceOwner()[fPopI];
                label patchi = mesh().boundaryMesh().whichPatch(fPopI);
                label lpfI = fPopI - mesh().boundaryMesh()[patchi].start();
                scalar sVf = phiPop.boundaryField()[patchi][lpfI]*
                    mesh().time().deltaTValue();
                sV[own] += sVf;
            }
        }

        forAll(popIndi, celli)
        {
            if (popIndi[celli])
            {
                oldVolgrow[celli] -= sV[celli];
            }
        }
    }

    forAll(mesh().cells(), celli)
    {
        if (popUpC[celli] == 1)
        {
            oldVolgrow[celli] = 0;
        }
    }

    growPopFields<scalar, fvPatchField, volMesh>
    (
        fPop,
        phiPop,
        popIndi,
        oldVolgrow
    );
    growPopFields<vector, fvPatchField, volMesh>
    (
        fPop,
        phiPop,
        popIndi,
        oldVolgrow
    );
    growPopFields<sphericalTensor, fvPatchField, volMesh>
    (
        fPop,
        phiPop,
        popIndi,
        oldVolgrow
    );
    growPopFields<symmTensor, fvPatchField, volMesh>
    (
        fPop,
        phiPop,
        popIndi,
        oldVolgrow
    );
    growPopFields<tensor, fvPatchField, volMesh>
    (
        fPop,
        phiPop,
        popIndi,
        oldVolgrow
    );
}


void Foam::fvMeshGIBChangersBase::correctV0()
{
    scalarField& V =
        const_cast<scalarField&>(dynamic_cast<const scalarField&>(mesh().V()));
    scalarField& V0 =
        const_cast<scalarField&>(dynamic_cast<const scalarField&>(mesh().V0()));
    surfaceScalarField& phiB = const_cast<surfaceScalarField&>(this->phiB());
    dimensionedScalar deltaT = mesh().time().deltaT();

    // Recompute V0
    V0 = V*(1 - scalarField(fvc::div(phiB))*deltaT.value());

    // Make negative V0 to small positive
    bool V0Negative = false;
    forAll(V0, celli)
    {
        if (V0[celli] < SMALL)
        {
            V0[celli] = SMALL;
            V0Negative = true;
        }
    }

    // The flag has to be synchronised otherwise in the equation solution for
    // correction the boundaries can't be synced and solver will freeze
    reduce(V0Negative, orOp<bool>());

    // Correct GIB fluxes to correspond to non-negative V0.
    // This happens for the badly snapped cells. Sometimes it can be avoided
    // by reducing the timestep.
    if (V0Negative)
    {
        // This function replaces modifyVolSqueezedCells and fixes the issue
        // with volumes not corresponding to fluxes. However, that does require
        // solving an equation (adds on computational time).
        negativeV0CorrectGIBFluxes(phiB, V, V0, deltaT);
    }
}


void Foam::fvMeshGIBChangersBase::negativeV0CorrectGIBFluxes
(
    surfaceScalarField& phiB,
    const scalarField& V,
    const scalarField& V0,
    const dimensionedScalar& deltaT
)
{
    // Initialize BCs list for GIBMeshPhiCorr to zero-gradient
    wordList GIBMeshPhiCorrTypes
    (
        phiB.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of GIBMeshPhiCorr to fixed-value for patches at which p is fixed
    forAll(phiB.boundaryField(), patchi)
    {
        if (phiB.boundaryField()[patchi].fixesValue())
        {
            GIBMeshPhiCorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField GIBMeshPhiCorr
    (
        IOobject
        (
            "GIBMeshPhiCorr",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensionedScalar
        (
            "GIBMeshPhiCorr",
            phiB.dimensions()*dimTime/dimLength,
            0.0
        ),
        GIBMeshPhiCorrTypes
    );
    mesh().schemes().setFluxRequired(GIBMeshPhiCorr.name());

    volScalarField divVolRequired
    (
        IOobject
        (
            "divVolRequired",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedScalar(dimless/dimTime, 0)
    );
    divVolRequired.primitiveFieldRef() = ((V - V0)/(V*deltaT.value()));

    fvScalarMatrix phiCorrEqn
    (
        fvm::laplacian(1.0/deltaT, GIBMeshPhiCorr)
     ==
        fvc::div(phiB) - divVolRequired
    );

    // Solver settings from cellDisplacement fvSolution sub-dictionary
    phiCorrEqn.solve(mesh().solution().solver("cellDisplacement"));
    phiB -= phiCorrEqn.flux();
}


void Foam::fvMeshGIBChangersBase::modifyVolSqueezedCells() const
{
    const scalarField& cV0 = mesh().V0();
    const scalarField& cV = mesh().V();

    DynamicList<label> dcV0(mesh().cells().size());
    DynamicList<label> dcV(mesh().cells().size());

    forAll(mesh().cells(), celli)
    {
        if (cV0[celli] < SMALL)
        {
            dcV0.append(celli);
        }
        if (cV[celli] < SMALL)
        {
            dcV.append(celli);
        }
    }
    dcV0.shrink();
    dcV.shrink();

    scalarField& V0 = const_cast<scalarField&>(cV0);
    scalarField& V = const_cast<scalarField&>(cV);

    const surfaceScalarField& cmPhi = mesh().phi();
    surfaceScalarField& mPhi = const_cast<surfaceScalarField&>(cmPhi);
    surfaceScalarField::Boundary& mPhibf = mPhi.boundaryFieldRef();

    forAll(dcV0, celli)
    {
        label gcI = dcV0[celli];

        if (V[gcI] < SMALL)
        {
            V[gcI] = SMALL;
        }
        V0[gcI] = V[gcI];

        const labelList& cFaces = mesh().cells()[gcI];
        forAll(cFaces, cFacei)
        {
            const label facei = cFaces[cFacei];
            if (facei < mesh().nInternalFaces())
            {
                mPhi[facei] = 0;
            }
            else
            {
                const label patchi = mesh().boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh().boundary()[patchi].patch();
                if (mPhibf[patchi].size() != 0)
                {
                    if (!isA<indirectPolyPatch>(pp))
                    {
                        const label pFacei = facei - pp.start();
                        mPhibf[patchi][pFacei] = 0;
                    }
                }
            }
        }
    }

    forAll(V, gcI)
    {
        if (V[gcI] < SMALL)
        {
            V[gcI] = SMALL;
            V0[gcI] = SMALL;
            const labelList& cFaces = mesh().cells()[gcI];

            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                if (facei < mesh().nInternalFaces())
                {
                    mPhi[facei] = 0;
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    const polyPatch& pp = mesh().boundary()[patchi].patch();
                    if (mPhibf[patchi].size() != 0)
                    {
                        if (!isA < indirectPolyPatch>(pp))
                        {
                            const label pFacei = facei - pp.start();
                            mPhibf[patchi][pFacei] = 0;
                        }
                    }
                }
            }
        }
    }

    if (debugMode_)
    {
        if (dcV0.size())
        {
            OFstream strV0
            (
                mesh().time().path()/
                "cells_V0_" + mesh().time().timeName() + ".obj"
            );
            meshTools::writeOBJ
            (
                strV0,
                mesh().cells(),
                mesh().faces(),
                mesh().points(),
                labelList(dcV0)
            );
        }

        if (dcV.size())
        {
            OFstream strV
            (
                mesh().time().path()/
                "cells_V_" + mesh().time().timeName() + ".obj"
            );
            meshTools::writeOBJ
            (
                strV,
                mesh().cells(),
                mesh().faces(),
                mesh().points(),
                labelList(dcV)
            );
        }
    }
}


void Foam::fvMeshGIBChangersBase::resetMeshFluxes()
{
    // Used to cancel the meshPhi with phi at boundary to make the relative 0,
    // but making that the mesh fluxes are inconsistent.
    const surfaceScalarField& phiBoundary = phiB();
    const surfaceScalarField& cmPhi = mesh().phi();
    surfaceScalarField& mPhi = const_cast<surfaceScalarField&>(cmPhi);
    mPhi = phiBoundary;
    volVectorField& U = mesh().lookupObjectRef<volVectorField>("U");

    U.correctBoundaryConditions();

    // Not consistent
    U.oldTime().correctBoundaryConditions();
}


void Foam::fvMeshGIBChangersBase::correctVelocityFlux()
{
    surfaceScalarField& phi = mesh().lookupObjectRef<surfaceScalarField>("phi");
    const volVectorField& U = mesh().lookupObject<volVectorField>("U");
    phi = fvc::interpolate(U) & mesh().Sf();
}


void Foam::fvMeshGIBChangersBase::visCells()
{
    masterFCells_ = 0;
    slaveFCells_ = 0;

    const labelList& fc1 = mesh().boundary()[masterGIB_].faceCells();
    const labelList& fc2 = mesh().boundary()[slaveGIB_].faceCells();
    forAll(fc1, fcI)
    {
        const label fc1I = fc1[fcI];
        if (fc1I < mesh().nInternalFaces())
        {
            masterFCells_[fc1I] = 1;
        }
        else
        {
            masterFCells_[fc1I] = 1;
        }
    }
    forAll(fc2, fcI)
    {
        const label fc2I = fc2[fcI];
        if (fc2I < mesh().nInternalFaces())
        {
            slaveFCells_[fc2I] = 1;
        }
        else
        {
            slaveFCells_[fc2I] = 1;
        }
    }
}


void Foam::fvMeshGIBChangersBase::faceCellsVisDebug()
{
    if (debugMode_)
    {
        visCells();
        writeScalarField("masterCells", masterFCells_);
        writeScalarField("slaveCells", slaveFCells_);
    }

    writeScalarField("cRegion", cRegion());
}


void Foam::fvMeshGIBChangersBase::findGIBPatches()
{
    patchGIBPtr_ = new boolList(mesh().boundary().size(), false);
    boolList& patchGIB = *patchGIBPtr_;
    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& poly = mesh().boundary()[patchi].patch();
        if (isA<indirectPolyPatch>(poly))
        {
            patchGIB[patchi] = true;

            const indirectPolyPatch& inPoly =
                refCast<const indirectPolyPatch>
                (
                    mesh().boundary()[patchi].patch()
                );

            const word& inPolyType = inPoly.indirectPolyPatchType();
            if (inPolyType == "master")
            {
                masterGIB_ = patchi;
            }
            else if (inPolyType == "slave")
            {
                slaveGIB_ = patchi;
            }
        }
    }
    if (masterGIB_ < 0 || slaveGIB_ < 0)
    {
        FatalErrorInFunction
            << "GIB master and slave patches were not found in mesh."
            << nl << exit(FatalError);
    }
}


void Foam::fvMeshGIBChangersBase::findGIBPatches(const word& zoneName)
{
    deleteDemandDrivenData(patchGIBPtr_);
    patchGIBPtr_ = new boolList(mesh().boundary().size(), false);
    boolList& patchGIB = *patchGIBPtr_;
    forAll(mesh().boundary(), patchi)
    {
        const polyPatch& poly = mesh().boundary()[patchi].patch();
        if (isA<indirectPolyPatch>(poly))
        {
            const indirectPolyPatch& inPoly =
                refCast<const indirectPolyPatch>
                (
                    mesh().boundary()[patchi].patch()
                );

            const label zoneId = inPoly.zoneId();
            if (mesh().faceZones()[zoneId].name() == zoneName)
            {
                patchGIB[patchi] = true;
                const word& inPolyType = inPoly.indirectPolyPatchType();
                if (inPolyType == "master")
                {
                    masterGIB_ = patchi;
                }
                else if (inPolyType == "slave")
                {
                    slaveGIB_ = patchi;
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::modifyRegionLabels(labelList& cellIndi) const
{
    const label fFluid = 0;
    const label fSolid = 1;
    labelList cellIndiBU = cellIndi;

    DynamicList<label> patchRegionLabel(mesh().boundary().size());

    forAll(region0Patch_, eI)
    {
        word substring = region0Patch_[eI];

        forAll(mesh().boundary(), patchi)
        {
            const fvPatch& cPatch(mesh().boundary()[patchi]);
            word patchName(cPatch.name());

            if (patchName.find(substring, 0) != string::npos)
            {
                const labelList& inCells = mesh().boundary()[patchi].faceCells();
                if (inCells.size()!=0)
                {
                    patchRegionLabel.append(cellIndiBU[inCells[0]]);
                }
            }
        }
    }
    patchRegionLabel.shrink();

    labelListList patchPro(Pstream::nProcs(), patchRegionLabel);

    Pstream::allGatherList(patchPro);

    forAll(cellIndi, celli)
    {
        bool found = false;
        forAll(patchPro, procI)
        {
            const labelList labels = patchPro[procI];
            forAll(labels, lI)
            {
                if (labels[lI] == cellIndiBU[celli])
                {
                    found = true;
                }
            }
        }
        if (!found)
        {
            cellIndi[celli] = fSolid;
        }
        else
        {
            cellIndi[celli] = fFluid;
        }
    }

    checkRegion(cellIndi);
}


bool Foam::fvMeshGIBChangersBase::isIdenticalInterface() const
{
    bool comparison = true;

    if (fl() != fl0())
    {
        comparison = false;
    }

    reduce(comparison, andOp<bool>());

    return comparison;
}


void Foam::fvMeshGIBChangersBase::smoothInternalBasePoints
(
    pointField& basePatchPoints,
    const labelList& pointAdd
) const
{
    const pointField& basePoints = *basePoints_;

    // Store the initial location of tge base interface in global addressing
    pointField gbasePatchPointsInit(basePoints);
    forAll(basePatchPoints, pointi)
    {
        gbasePatchPointsInit[pointAdd[pointi]] = basePatchPoints[pointi];
    }
    syncPoints(gbasePatchPointsInit);

    // Mark points which are on faces connected on the boundary
    const labelList& interFaces = fl();
    const boolList& bPoints = boundaryPoints();

    // Smooth base patch with relaxation
    for (int i = 1; i < 4; i++)
    {
        pointField gbasePatchPoints(basePoints);
        forAll(basePatchPoints, pointi)
        {
            gbasePatchPoints[pointAdd[pointi]] = basePatchPoints[pointi];
        }

        // Construct patch
        indirectPrimitivePatch smoothPatch
        (
            IndirectList<face>(mesh().faces(), interFaces),
            gbasePatchPoints
        );

        PrimitivePatchInterpolation<indirectPrimitivePatch> pInterC
        (
            smoothPatch
        );

        tmp<pointField> smoothFaceAveraget =
            pInterC.faceToPointInterpolate(smoothPatch.faceCentres());
        pointField& smoothFaceAverage = smoothFaceAveraget.ref();

        const labelList& globalAddressing = smoothPatch.meshPoints();
        forAll(smoothFaceAverage, pointi)
        {
            const label gpI = globalAddressing[pointi];
            if (bPoints[gpI])
            {
                smoothFaceAverage[pointi] = gbasePatchPoints[gpI];
            }
        }

        pointField gsmoothFaceAverage(basePoints);
        forAll(smoothFaceAverage, pointi)
        {
            gsmoothFaceAverage[globalAddressing[pointi]] =
                smoothFaceAverage[pointi];
        }
        syncPoints(gsmoothFaceAverage);

        correctConstraintPatches(gsmoothFaceAverage);

        // Update the baseMesh points using relaxation
        forAll(basePatchPoints, pointi)
        {
            const label gpI = pointAdd[pointi];
            basePatchPoints[pointi] =
                0.5*gbasePatchPointsInit[gpI]
              + 0.5*gsmoothFaceAverage[gpI];
        }
    }
}


void Foam::fvMeshGIBChangersBase::updateSolidCellZone() const
{
    const labelList& cReg = cRegion();
    const cellZone& cZone = mesh().cellZones()[findOrCreateCellZone()];
    cellZone& ccZone = const_cast<cellZone&>(cZone);

    DynamicList<label> zoneCL(cReg.size());

    forAll(cReg, celli)
    {
        if (cReg[celli] != 0)
        {
            zoneCL.append(celli);
        }
    }
    zoneCL.shrink();
    ccZone.clearAddressing();
    const labelList c(zoneCL);
    ccZone = c;
}


void Foam::fvMeshGIBChangersBase::storeOldTimes()
{
    StoreOldTimeFields<scalar, fvPatchField, volMesh>(mesh());
    StoreOldTimeFields<vector, fvPatchField, volMesh>(mesh());
    StoreOldTimeFields<sphericalTensor, fvPatchField, volMesh>(mesh());
    StoreOldTimeFields<symmTensor, fvPatchField, volMesh>(mesh());
    StoreOldTimeFields<tensor, fvPatchField, volMesh>(mesh());
    StoreOldTimeFields<scalar, fvsPatchField, surfaceMesh>(mesh());
    StoreOldTimeFields<vector, fvsPatchField, surfaceMesh>(mesh());
    StoreOldTimeFields<sphericalTensor, fvsPatchField, surfaceMesh>(mesh());
    StoreOldTimeFields<symmTensor, fvsPatchField, surfaceMesh>(mesh());
    StoreOldTimeFields<tensor, fvsPatchField, surfaceMesh>(mesh());
}


void Foam::fvMeshGIBChangersBase::correctBCs()
{
    // Only sync parallel BCs as there are issues correcting some BCs
    // at an intermediate stage
    CorrectProcessorBCs<scalar>(mesh());
    CorrectProcessorBCs<vector>(mesh());
    CorrectProcessorBCs<sphericalTensor>(mesh());
    CorrectProcessorBCs<symmTensor>(mesh());
    CorrectProcessorBCs<tensor>(mesh());

    if (correctBCs_)
    {
        if (mesh().foundObject<volScalarField>("k"))
        {
            volScalarField& k = mesh().lookupObjectRef<volScalarField>("k");
            forAll(k, celli)
            {
                if (k[celli] < SMALL)
                {
                    k[celli] = SMALL;
                }
            }
            volScalarField::Boundary& kbf = k.boundaryFieldRef();
            forAll(kbf, patchi)
            {
                forAll(kbf[patchi], pFacei)
                {
                    if (kbf[patchi][pFacei] < SMALL)
                    {
                        kbf[patchi][pFacei] = SMALL;
                    }
                }
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::explicitLimitingFields()
{
    if (limitNegFieldNames_.size())
    {
        forAll(limitNegFieldNames_, fieldI)
        {
            word substring = limitNegFieldNames_[fieldI];

            if (mesh().foundObject<volScalarField>(substring))
            {
                volScalarField& field =
                    mesh().lookupObjectRef<volScalarField>(substring);

                label nT = field.nOldTimes();

                forAll(field, celli)
                {
                    if (field[celli] < 0)
                    {
                        field[celli] = 0;
                    }
                }
                if (nT != 0)
                {
                    forAll(field.oldTime(), celli)
                    {
                        if (field.oldTime()[celli] < 0)
                        {
                            field.oldTime()[celli] = 0;
                        }
                    }
                }
            }
            else
            {
                FatalErrorInFunction << abort(FatalError);
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::clearOutGIBData()
{
    deleteDemandDrivenData(multInterFacesPtr_);
    deleteDemandDrivenData(flPtr_);
    deleteDemandDrivenData(fl0Ptr_);
    deleteDemandDrivenData(fmPtr_);
    deleteDemandDrivenData(fm0Ptr_);
    deleteDemandDrivenData(interPointsPtr_);
    deleteDemandDrivenData(interPoints0Ptr_);
    deleteDemandDrivenData(phiBPtr_);
    deleteDemandDrivenData(phiTrPtr_);
    deleteDemandDrivenData(phiTrSPtr_);
    deleteDemandDrivenData(phiTrGPtr_);
    deleteDemandDrivenData(faceIndicatorPtr_);
    deleteDemandDrivenData(faceIndicator0Ptr_);
    deleteDemandDrivenData(cRegionPtr_);
    deleteDemandDrivenData(cRegion0Ptr_);
    deleteDemandDrivenData(recPoints0Ptr_);
    deleteDemandDrivenData(recAllPoints0Ptr_);
    deleteDemandDrivenData(hitPointPtr_);
    deleteDemandDrivenData(closeBoundaryPointsPtr_);
    deleteDemandDrivenData(markedBoundaryPointsPtr_);
    deleteDemandDrivenData(normal2DEdgesPtr_);
    deleteDemandDrivenData(boundaryPointsRevPtr_);
    deleteDemandDrivenData(fullSnapCellsPtr_);
    deleteDemandDrivenData(fullSnapCellPointsPtr_);
    deleteDemandDrivenData(popCellPointsPtr_);
    deleteDemandDrivenData(popSPointsPtr_);
    deleteDemandDrivenData(popGPointsPtr_);
    deleteDemandDrivenData(popUpCellsPtr_);
    deleteDemandDrivenData(meshQualityPointsPtr_);
    deleteDemandDrivenData(patchGIBPtr_);
    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(concavePointsPtr_);
    deleteDemandDrivenData(concaveEdgesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangersBase::fvMeshGIBChangersBase
(
    fvMesh& mesh,
    const word& typeN
)
:
    fvMeshGIBChanger(mesh),
    dynamicMeshCoeffs_(dict().subDict(typeN + "Coeffs")),
    basePoints_(nullptr),
    baseCC_(nullptr),
    baseSf_(nullptr),
    baseCf_(nullptr),
    triName_(dynamicMeshCoeffs_.lookup("triSurfaceName")),
    patchGIBPtr_(nullptr),
    masterGIB_(-1),
    slaveGIB_(-1),
    cRegionPtr_(nullptr),
    cRegion0Ptr_(nullptr),
    closeBoundaryPointsPtr_(nullptr),
    markedBoundaryPointsPtr_(nullptr),
    normal2DEdgesPtr_(nullptr),
    boundaryPointsRevPtr_(nullptr),
    ibMeshPtr_(nullptr),
    flPtr_(nullptr),
    fmPtr_(nullptr),
    interPointsPtr_(nullptr),
    interPoints0Ptr_(nullptr),
    boundaryPointsPtr_(nullptr),
    multInterFacesPtr_(nullptr),
    concavePointsPtr_(nullptr),
    concaveEdgesPtr_(nullptr),
    fullSnapCellsPtr_(nullptr),
    fullSnapCellPointsPtr_(nullptr),
    popCellPointsPtr_(nullptr),
    popSPointsPtr_(nullptr),
    popGPointsPtr_(nullptr),
    hitIndexPtr_(nullptr),
    hitPointPtr_(nullptr),
    popUpCellsPtr_(nullptr),
    faceIndicatorPtr_(nullptr),
    faceIndicator0Ptr_(nullptr),
    masterFCells_(mesh.cells().size(), false),
    slaveFCells_(mesh.cells().size(), false),
    region0Patch_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("region0Patch", wordList(1,"inlet"))
    ),
    limitNegFieldNames_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("limitNegFieldNames", wordList())
    ),
    popSubSteps_
    (
        dynamicMeshCoeffs_.lookupOrDefault<label>
        ("popSubSteps", 10)
    ),
    maxNonOrth_
    (
        cos
        (
            degToRad
            (
                dynamicMeshCoeffs_.lookupOrDefault<scalar>
                ("maxNonOrthogonality", 80)
            )
        )
    ),
    unsnapVar_
    (
        dynamicMeshCoeffs_.lookupOrDefault<scalar>
        ("unsnapCellsPercentage", 0.3)
    ),
    fl0Ptr_(nullptr),
    fm0Ptr_(nullptr),
    phiBPtr_(nullptr),
    phiTrPtr_(nullptr),
    phiTrSPtr_(nullptr),
    phiTrGPtr_(nullptr),
    debugMode_(dynamicMeshCoeffs_.lookupOrDefault<bool>("debug", false)),
    correctBCs_(dynamicMeshCoeffs_.lookupOrDefault<bool>("correctBCs", true)),
    pyrPrismFlip_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("flipPyramidsAndPrisms", false)
    ),
    boundPopValues_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("boundPopValues", true)
    ),
    checkInterface_(false),
    oldPoints_(mesh.points()),
    prevPoints_(mesh.points()),
    recPoints0Ptr_(nullptr),
    recAllPoints0Ptr_(nullptr),
    meshQualityPointsPtr_(nullptr)
{
    initialize();
}


Foam::fvMeshGIBChangersBase::fvMeshGIBChangersBase
(
    fvMesh& mesh,
    const dictionary dict
)
:
    fvMeshGIBChanger(mesh),
    dynamicMeshCoeffs_
    (
        dict.optionalSubDict(dict.lookup<word>("type") + "Coeffs")
    ),
    basePoints_(nullptr),
    baseCC_(nullptr),
    baseSf_(nullptr),
    baseCf_(nullptr),
    triName_(dynamicMeshCoeffs_.lookup("triSurfaceName")),
    patchGIBPtr_(nullptr),
    masterGIB_(-1),
    slaveGIB_(-1),
    cRegionPtr_(nullptr),
    cRegion0Ptr_(nullptr),
    closeBoundaryPointsPtr_(nullptr),
    markedBoundaryPointsPtr_(nullptr),
    normal2DEdgesPtr_(nullptr),
    boundaryPointsRevPtr_(nullptr),
    ibMeshPtr_(nullptr),
    flPtr_(nullptr),
    fmPtr_(nullptr),
    interPointsPtr_(nullptr),
    interPoints0Ptr_(nullptr),
    boundaryPointsPtr_(nullptr),
    multInterFacesPtr_(nullptr),
    concavePointsPtr_(nullptr),
    concaveEdgesPtr_(nullptr),
    fullSnapCellsPtr_(nullptr),
    fullSnapCellPointsPtr_(nullptr),
    popCellPointsPtr_(nullptr),
    popSPointsPtr_(nullptr),
    popGPointsPtr_(nullptr),
    hitIndexPtr_(nullptr),
    hitPointPtr_(nullptr),
    popUpCellsPtr_(nullptr),
    faceIndicatorPtr_(nullptr),
    faceIndicator0Ptr_(nullptr),
    masterFCells_(mesh.cells().size(), false),
    slaveFCells_(mesh.cells().size(), false),
    region0Patch_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("region0Patch", wordList(1,"inlet"))
    ),
    limitNegFieldNames_
    (
        dynamicMeshCoeffs_.lookupOrDefault<wordList>
        ("limitNegFieldNames", wordList())
    ),
    popSubSteps_
    (
        dynamicMeshCoeffs_.lookupOrDefault<label>
        ("popSubSteps", 10)
    ),
    maxNonOrth_
    (
        cos
        (
            degToRad
            (
                dynamicMeshCoeffs_.lookupOrDefault<scalar>
                ("maxNonOrthogonality", 80)
            )
        )
    ),
    unsnapVar_
    (
        dynamicMeshCoeffs_.lookupOrDefault<scalar>
        ("unsnapCellsPercentage", 0.3)
    ),
    fl0Ptr_(nullptr),
    fm0Ptr_(nullptr),
    phiBPtr_(nullptr),
    phiTrPtr_(nullptr),
    phiTrSPtr_(nullptr),
    phiTrGPtr_(nullptr),
    debugMode_(dynamicMeshCoeffs_.lookupOrDefault<bool>("debug", false)),
    correctBCs_(dynamicMeshCoeffs_.lookupOrDefault<bool>("correctBCs", true)),
    pyrPrismFlip_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("flipPyramidsAndPrisms", false)
    ),
    boundPopValues_
    (
        dynamicMeshCoeffs_.lookupOrDefault<bool>("boundPopValues", true)
    ),
    checkInterface_(false),
    oldPoints_(mesh.points()),
    prevPoints_(mesh.points()),
    recPoints0Ptr_(nullptr),
    recAllPoints0Ptr_(nullptr),
    meshQualityPointsPtr_(nullptr)
{
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangersBase::~fvMeshGIBChangersBase()
{
    deleteDemandDrivenData(basePoints_);
    deleteDemandDrivenData(baseCC_);
    deleteDemandDrivenData(baseSf_);
    deleteDemandDrivenData(baseCf_);
    deleteDemandDrivenData(ibMeshPtr_);
    deleteDemandDrivenData(hitIndexPtr_);

    clearOutGIBData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshGIBChangersBase::updateInit(const word& zoneName)
{
    findGIBPatches(zoneName);
    clearOutGIBData();

    const label zoneId = gibFaceZone();
    const faceZone& cfZone = mesh().faceZones()[zoneId];
    faceZone& fZone = const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(0);
    fm0Ptr_ = new boolList(0);

    tmp<pointField> snapPt = findSnappedPoints(false);
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList& fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    if (tangentialBoundaryMotion())
    {
        doTangentialBoundaryMotion(snapP);
    }
    correctConstraintPatches(snapP);

    mesh().moveGIBPoints(snapP);
    mesh().updateGIB();

    word cRegionName = "cRegion";
    if (zoneName != word::null)
    {
        cRegionName += "_"+zoneName;
    }

    writeScalarField(cRegionName, cRegion(), true);
}


void Foam::fvMeshGIBChangersBase::doUpdate(GIBMapping& mapCl, bool fromBase)
{
    const label zoneId = gibFaceZone();
    const faceZone& cfZone = mesh().faceZones()[zoneId];
    faceZone& fZone = const_cast<faceZone&>(cfZone);
    fl0Ptr_ = new labelList(fZone);
    fm0Ptr_ = new boolList(fZone.flipMap());

    tmp<pointField> snapPt = findSnappedPoints(fromBase);
    pointField& snapP = snapPt.ref();

    const labelList& fzAdd = fl();
    const boolList& fzFm = fm();
    fZone.resize(fzAdd.size());
    fZone.resetAddressing(fzAdd,fzFm);

    if (tangentialBoundaryMotion())
    {
        doTangentialBoundaryMotion(snapP);
    }

    correctConstraintPatches(snapP);

    mesh().moveGIBPoints(snapP);

    checkInterface_ = isIdenticalInterface();

    if (!checkInterface_)
    {
        mesh().updateGIB();

        mapCl.mapBcs();

        popShrinkFields();

        correctV0();

        popGrowFields();
    }
    else
    {
        mesh().storeGIBFields();
        correctV0();
    }

    oldPoints_ = recAllPoints0();

    resetMeshFluxes();

    writeProblematicCells();

    faceCellsVisDebug();

    correctBCs();

    explicitLimitingFields();
}


void Foam::fvMeshGIBChangersBase::writeGeometry(const fileName& surfaceName)
{
    const fvPatch& mPatch = mesh().boundary()[masterId()];

    MeshedSurface<face> surface(mPatch.patch(), true);

    if (Pstream::master())
    {

        fileName globalCasePath
        (
            surfaceName.isAbsolute()
          ? surfaceName
          : (
                mesh().time().processorCase()
              ? mesh().time().rootPath()
               /mesh().time().globalCaseName()/surfaceName
              : mesh().time().path()/surfaceName
            )
        );
        globalCasePath.clean();

        if (!exists(globalCasePath.path()))
        {
            Info<< "Path to surface does not exist." << endl;
            Info<< "Creating path: "
                 << globalCasePath.path()
                 << endl;

            mkDir(globalCasePath.path());
        }

        surface.write(globalCasePath);
    }
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::fl() const
{
    if (!flPtr_)
    {
        makeFl();
    }

    return *flPtr_;
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::fl0() const
{
    return *fl0Ptr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::fm() const
{
    if (!fmPtr_)
    {
        makeFlipMap();
    }

    return *fmPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::fm0() const
{
    return *fm0Ptr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::interP() const
{
    if (!interPointsPtr_)
    {
        makeInterPoints();
    }

    return *interPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::interP0() const
{
    if (!interPoints0Ptr_)
    {
        makeInterPoints0();
    }
    return *interPoints0Ptr_;
}


const Foam::surfaceScalarField& Foam::fvMeshGIBChangersBase::phiB() const
{
    if (!phiBPtr_)
    {
        makePhiB();
    }

    return *phiBPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshGIBChangersBase::phiTr() const
{
    if (!phiTrPtr_)
    {
        makePhiTr();
    }

    return *phiTrPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshGIBChangersBase::phiTrS() const
{
    if (!phiTrSPtr_)
    {
        makePhiTrS();
    }

    return *phiTrSPtr_;
}


const Foam::surfaceScalarField& Foam::fvMeshGIBChangersBase::phiTrG() const
{
    if (!phiTrGPtr_)
    {
        makePhiTrG();
    }

    return *phiTrGPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::faceIndicator() const
{
    if (!faceIndicatorPtr_)
    {
        makeFaceIndicator();
    }

    return *faceIndicatorPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::faceIndicator0() const
{
    if (!faceIndicator0Ptr_)
    {
        makeFaceIndicator0();
    }

    return *faceIndicator0Ptr_;
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::cRegion() const
{
    if (!cRegionPtr_)
    {
        makecRegion();
    }

    return *cRegionPtr_;
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::cRegion0() const
{
    if (!cRegion0Ptr_)
    {
        makecRegion0();
    }

    return *cRegion0Ptr_;
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::fullSnapCells() const
{
    if (!fullSnapCellsPtr_)
    {
        makeFullSnapCells();
    }

    return *fullSnapCellsPtr_;
}


const Foam::labelList& Foam::fvMeshGIBChangersBase::fullSnapCellPoints() const
{
    if (!fullSnapCellPointsPtr_)
    {
        makeFullSnapCellPoints();
    }

    return *fullSnapCellPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::popCellPoints() const
{
    if (!popCellPointsPtr_)
    {
        makePopCellPoints();
    }

    return *popCellPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::popSPoints() const
{
    if (!popSPointsPtr_)
    {
        makePopSPoints();
    }

    return *popSPointsPtr_;
}

const Foam::boolList& Foam::fvMeshGIBChangersBase::popGPoints() const
{
    if (!popGPointsPtr_)
    {
        makePopGPoints();
    }

    return *popGPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::concavePoints() const
{
    if (!concavePointsPtr_)
    {
        makeConcavePoints();
    }

    return *concavePointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::concaveEdges() const
{
    if (!concaveEdgesPtr_)
    {
        makeConcavePoints();
    }

    return *concaveEdgesPtr_;
}


const Foam::labelListList&
Foam::fvMeshGIBChangersBase::closeBoundaryPoints() const
{
    if (!closeBoundaryPointsPtr_)
    {
        makeCloseBoundaryPoints();
    }

    return *closeBoundaryPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::markedBoundaryPoints() const
{
    if (!markedBoundaryPointsPtr_)
    {
        makeMarkedBoundaryPoints();
    }

    return *markedBoundaryPointsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::boundaryPointsRev() const
{
    if (!boundaryPointsRevPtr_)
    {
        makeBoundaryPointsRev();
    }

    return *boundaryPointsRevPtr_;
}


const Foam::pointField& Foam::fvMeshGIBChangersBase::recPoints0() const
{
    if (!recPoints0Ptr_)
    {
        makeRecPoints0();
    }

    return *recPoints0Ptr_;
}


const Foam::pointField& Foam::fvMeshGIBChangersBase::recAllPoints0() const
{
    if (!recAllPoints0Ptr_)
    {
        makeRecAllPoints0();
    }

    return *recAllPoints0Ptr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::popUpCells() const
{
    if (!popUpCellsPtr_)
    {
        makePopUpCells();
    }

    return *popUpCellsPtr_;
}


const Foam::boolList& Foam::fvMeshGIBChangersBase::boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        makeBoundaryPoints();
    }

    return *boundaryPointsPtr_;
}


// ************************************************************************* //
