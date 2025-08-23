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
    (c) 2022-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "equationRelax.H"
#include "cellQuality/cellQuality.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

void Foam::equationRelax::initialize()
{
    if (isTransient() && localRelax_)
    {
        //- Safety check error
        FatalErrorInFunction
            << "Local relaxation is not supported for transient runs"
            << abort(FatalError);
    }
    if (localRelax_)
    {
        const dictionary& locRelaxDict = dict_.subDict("localRelaxCoeffs");
        relaxType_ = locRelaxDict.lookupOrDefault<word>("relaxType", "Courant");
        if (relaxType_ == "Courant")
        {
            Info<< "Steady-state run ";
            Info<< "using Courant-based URF relaxation" << endl;
            if (locRelaxDict.found("minCourant"))
            {
                minC_ = Function1<scalar>::New("minCourant", locRelaxDict);
            }
            else
            {
                minC_.reset
                (
                    new Function1Types::Constant<scalar>("minCourant", 0.1)
                );
            }
            if (locRelaxDict.found("maxCourant"))
            {
                maxC_ = Function1<scalar>::New("maxCourant", locRelaxDict);
            }
            else
            {
                maxC_.reset
                (
                    new Function1Types::Constant<scalar>("maxCourant", 1e+5)
                );
            }
            if (locRelaxDict.found("minCourantRelax"))
            {
                minCRelax_ =
                    Function1<scalar>::New("minCourantRelax", locRelaxDict);
            }
            else
            {
                minCRelax_.reset
                (
                    new Function1Types::Constant<scalar>
                    (
                         "minCourantRelax", 0.999999
                    )
                );
            }
            if (locRelaxDict.found("maxCourantRelax"))
            {
                maxCRelax_ =
                    Function1<scalar>::New("maxCourantRelax", locRelaxDict);
            }
            else
            {
                maxCRelax_.reset
                (
                    new Function1Types::Constant<scalar>
                    (
                        "maxCourantRelax", 0.95
                    )
                );
            }
            calculateCellDeltas();
        }
        else
        {
            FatalErrorInFunction
                << "Unknown relaxation type " << relaxType_ << endl
                << exit(FatalError);
        }
    }
    if (meshQualityRelax_ < 1)
    {
        const dictionary& meshQdict =
            dict_.subDict("meshQualityControls");
        nonOrthoThreshold_ =
            meshQdict.lookupOrDefault<scalar>("nonOrthoThreshold", 60.0);
        skewnessThreshold_ =
            meshQdict.lookupOrDefault<scalar>("skewnessThreshold", 0.9);

        // get the mesh quality field
        badCellQualityMarker_ = badCellQualityMarker();
    }
    if (localCellProcRelax_ < 1)
    {
        badDecompositionCellMarker_ = badDecompositionCellMarker();
    }
    if (localCellProcAMIRelax_ < 1)
    {
        procAndAMICellMarker_ = cellsWithProcAndAMIFaces();
    }
}


void Foam::equationRelax::nonUniformMomentumRelax
(
    fvBlockMatrix<vector>& bUEq
) const
{
    tmp<vectorField> locR = computeLocalRelaxFactor(bUEq.psi());
    bUEq.relax(locR());
}


tmp<Foam::vectorField> Foam::equationRelax::computeLocalRelaxFactor
(
    const volVectorField& U
) const
{
    tmp<vectorField> trelaxFactor
    (
        new vectorField
        (
            mesh_.nCells(),
            vector::one
        )
    );
    if (relaxType_ == "Courant")
    {
        scalar minC = minC_->value(mesh_.time().timeIndex());
        scalar maxC = maxC_->value(mesh_.time().timeIndex());
        scalar minCourantRelax = minCRelax_->value(mesh_.time().timeIndex());
        scalar maxCourantRelax = maxCRelax_->value(mesh_.time().timeIndex());
        tmp<vectorField> tvCo = computeCourantNo(U);
        const vectorField& vCo = tvCo.ref();
        vectorField& relaxFactor = trelaxFactor.ref();
        // Compute the relax depending in each direction
        forAll(relaxFactor, cI)
        {
            const vector& vCoi = vCo[cI];
            forAll(vCoi, dI)
            {
                const scalar c = vCoi[dI];
                scalar relax(0.0);
                if (c>maxC)
                {
                    relax = maxCourantRelax;
                }
                else if (c<minC)
                {
                    relax = minCourantRelax;
                }
                else
                {
                    relax = maxCourantRelax +
                        (minCourantRelax-maxCourantRelax)*(maxC-c)/(maxC-minC);
                }
                relaxFactor[cI][dI] = relax;
            }
        }
    }
    if (false && mesh_.time().outputTime())
    {
        volVectorField vRelVis
        (
            IOobject
            (
                "vRelVis",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector(dimless, vector(0, 0, 0)),
            zeroGradientFvPatchScalarField::typeName
        );
        vRelVis.primitiveFieldRef() = trelaxFactor();
        vRelVis.correctBoundaryConditions();
        vRelVis.write();
    }
    return trelaxFactor;
}


tmp<Foam::vectorField> Foam::equationRelax::computeCourantNo
(
    const volVectorField& U
) const
{
    tmp<vectorField> tCoi
    (
        new vectorField
        (
            mesh_.nCells(),
            vector::one
        )
    );
    vectorField& Coi = tCoi.ref();

    //- Sanity
    if (!pcdelta_.valid())
    {
            FatalErrorInFunction
        << "Cell deltas are not computed"
        << exit(FatalError);
    }

    //- calculate the max deltas for xyz
    const vectorField& pcdelta = pcdelta_();
    forAll(Coi, celli)
    {
        const vector& dx = pcdelta[celli];
        Coi[celli].x() = mag(U[celli].x())/dx.x();
        Coi[celli].y() = mag(U[celli].y())/dx.y();
        Coi[celli].z() = mag(U[celli].z())/dx.z();
    }

    if (false && mesh_.time().outputTime())
    {
        volVectorField vCoVis
        (
            IOobject
            (
                "vCoVis",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector(dimless, vector(0, 0, 0)),
            zeroGradientFvPatchScalarField::typeName
        );
        vCoVis.primitiveFieldRef() = Coi;
        vCoVis.correctBoundaryConditions();
        vCoVis.write();
    }
    return tCoi;
}


void Foam::equationRelax::calculateCellDeltas()
{
    pcdelta_.reset(new vectorField(mesh_.nCells(), Zero));
    pcdelta_.reset(new vectorField(mesh_.nCells(), Zero));
    vectorField minC(mesh_.nCells(), vector::max);
    vectorField maxC(mesh_.nCells(), vector::min);
    const pointField& points = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();

    //- calculate min-max coordinate for each cell
    forAll(points, pointi)
    {
        const vector& p = points[pointi];
        const labelList& pointCellsI = pointCells[pointi];

        forAll(pointCellsI, pcI)
        {
            const label celli = pointCellsI[pcI];
            minC[celli] = min(minC[celli], p);
            maxC[celli] = max(maxC[celli], p);
        }
    }
    vectorField& pcdelta = pcdelta_();
    forAll(pcdelta, celli)
    {
        pcdelta[celli] = maxC[celli] - minC[celli];
    }
}


Foam::labelList Foam::equationRelax::badCellQualityMarker(bool printInfo) const
{
    DynamicList<label> badQualityCellMap;
    if (meshQualityRelax_<1.0)
    {
        cellQuality cq(mesh_);
        tmp<scalarField> tcNonOrtho = cq.nonOrthogonality();
        tmp<scalarField> tcSkew = cq.skewness();
        const scalarField cNonOrtho = tcNonOrtho();
        const scalarField cSkew = tcSkew();

        forAll(cNonOrtho, cI)
        {
            if
            (
                (cNonOrtho[cI]>nonOrthoThreshold_)
              ||(cSkew[cI] >skewnessThreshold_)
            )
            {
                badQualityCellMap.append(cI);
            }
        }
        badQualityCellMap.shrink();

        if (printInfo)
        {
            label nCells = badQualityCellMap.size();
            reduce(nCells, sumOp<label>());

            Info<< "Cells to be relaxed due to quality: " << nCells
                 << nl << endl;
        }
    }
    return labelList(badQualityCellMap);
}


Foam::labelList Foam::equationRelax::badDecompositionCellMarker
(
    bool printInfo
) const
{
    DynamicList<label> badDecompositionCellMap;
    if (localCellProcRelax_ < 1.0)
    {
        label nSol = mesh_.nSolutionD();

        if (nSol == 3)
        {
            const cellList& cells = mesh_.cells();
            forAll(cells, cI)
            {
                label nInternal = 0;
                bool wallCell = false;
                forAll(cells[cI], fI)
                {
                    label gfI = cells[cI][fI];
                    if (gfI < mesh_.nInternalFaces())
                    {
                        nInternal++;
                    }
                    else
                    {
                        label patchi = mesh_.boundaryMesh().whichPatch(gfI);
                        if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
                        {
                            wallCell = true;
                        }
                    }
                }
                if
                (
                    (nInternal<=2) && wallCell
                )
                {
                    badDecompositionCellMap.append(cI);
                }
            }
            badDecompositionCellMap.shrink();

            if (printInfo)
            {
                label cellMapSize(badDecompositionCellMap.size());
                reduce(cellMapSize, sumOp<label>());

                Info<< "Wall cells with only one or two internal faces: "
                     << cellMapSize << nl << endl;
            }
        }
    }
    return labelList(badDecompositionCellMap);
}


Foam::labelList Foam::equationRelax::cellsWithProcAndAMIFaces
(
    bool printInfo
) const
{
    DynamicList<label> cellMap;

    if (localCellProcAMIRelax_<1.0)
    {
        //first mark all the procCells
        boolList isProcCell(mesh_.nCells(), false);
        boolList isProcAMICell(mesh_.nCells(), false);

        const cellList& cells = mesh_.cells();
        forAll(mesh_.boundary(), pI)
        {
            if (isA<processorFvPatch>(mesh_.boundary()[pI]))
            {
                const labelList& fcs = mesh_.boundary()[pI].faceCells();
                forAll(fcs, cI)
                {
                    const label gcI = fcs[cI];
                    if (!isProcCell[gcI])
                    {
                        isProcCell[gcI] = true;
                        if (!isProcAMICell[gcI])
                        {
                            forAll(cells[gcI], fI)
                            {
                                label gfI = cells[gcI][fI];
                                if (!(gfI < mesh_.nInternalFaces()))
                                {
                                    const label pII =
                                        mesh_.boundaryMesh().whichPatch(gfI);
                                    if
                                    (
                                        isA<cyclicAMIFvPatch>
                                        (
                                            mesh_.boundary()[pII]
                                        )
                                    )
                                    {
                                        isProcAMICell[gcI] = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        forAll(isProcAMICell, cI)
        {
            if (isProcAMICell[cI])
            {
                cellMap.append(cI);
            }
        }
        if (false)
        {
            volScalarField visProcAMIcells
            (
                IOobject
                (
                    "visProcAMIcells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 0),
                "zeroGradient"
            );
            forAll(visProcAMIcells, cI)
            {
                if (isProcAMICell[cI]) visProcAMIcells[cI] = 1;
            }
            visProcAMIcells.write();
        }
        cellMap.shrink();

        if (printInfo)
        {
            label cellMapSize(cellMap.size());
            reduce(cellMapSize, sumOp<label>());

            Info<< "Cells with both AMI and processor faces: "
                 << cellMapSize << nl << endl;
        }
    }
    return labelList(cellMap, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationRelax::equationRelax
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool fromSolverObject
)
:
    mesh_(mesh),
    dict_(dict),
    timeMode_
    (
        fromSolverObject?
        dict.lookupOrDefault<word>("timeMode", "physical") : "pseudoTime"
    ),
    localRelax_(dict.lookupOrDefault<Switch>("localRelax", false)),
    relaxType_(word::null),
    minC_(),
    maxC_(),
    minCRelax_(),
    maxCRelax_(),
    procRelax_(dict.lookupOrDefault<scalar>("procRelax", 1.0)),
    wallRelax_(dict.lookupOrDefault<scalar>("wallRelax", 1.0)),
    meshQualityRelax_(dict.lookupOrDefault<scalar>("meshQualityRelax", 1.0)),
    nonOrthoThreshold_(0.0),
    skewnessThreshold_(0.0),
    badCellQualityMarker_(0),
    localCellProcRelax_
    (
        dict.lookupOrDefault<scalar>("localCellProcRelax", 1.0)
    ),
    badDecompositionCellMarker_(0),
    localCellProcAMIRelax_
    (
        dict.lookupOrDefault<scalar>("localCellProcAMIRelax", 1.0)
    ),
    procAndAMICellMarker_(0),
    pcdelta_(nullptr)
{
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationRelax::momentumRelax(fvBlockMatrix<vector>& bUEq) const
{
    bUEq.boundaryRelax();
    if (procRelax_ < 1)
    {
        boolList isProcBoundary(mesh_.boundary().size(), false);
        forAll(isProcBoundary, pI)
        {
            const fvPatch& fvp =  mesh_.boundary()[pI];
            if (isA<processorFvPatch>(fvp)) isProcBoundary[pI] = true;
        }
        bUEq.boundaryRelax(procRelax_, isProcBoundary);
    }
    if (wallRelax_ < 1)
    {
        boolList isWallBoundary(mesh_.boundary().size(), false);
        forAll(isWallBoundary, pI)
        {
            const fvPatch& fvp =  mesh_.boundary()[pI];
            if (isA<wallFvPatch>(fvp)) isWallBoundary[pI] = true;
        }
        bUEq.boundaryRelax(wallRelax_, isWallBoundary);
    }
    if (meshQualityRelax_ < 1)
    {
        bUEq.relax(meshQualityRelax_, badCellQualityMarker_);
    }
    if (localCellProcRelax_ < 1)
    {
        bUEq.relax(localCellProcRelax_, badDecompositionCellMarker_);
    }
    if (localCellProcAMIRelax_ < 1)
    {
        bUEq.relax(localCellProcAMIRelax_, procAndAMICellMarker_);
    }
    bUEq.relax();
    if (localRelax_)
    {
        nonUniformMomentumRelax(bUEq);
    }
}


void Foam::equationRelax::updateCellMarkers()
{
    if (localRelax_)
    {
        if (relaxType_ == "Courant") calculateCellDeltas();
    }
    if (meshQualityRelax_ < 1)
    {
        badCellQualityMarker_ = badCellQualityMarker(false);
    }
    if (localCellProcRelax_ < 1)
    {
        badDecompositionCellMarker_ = badDecompositionCellMarker(false);
    }
    if (localCellProcAMIRelax_ < 1)
    {
        procAndAMICellMarker_ = cellsWithProcAndAMIFaces(false);
    }
}

// ************************************************************************* //
