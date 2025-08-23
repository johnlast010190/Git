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

#include "flowSolver.H"
#include "solverObjects/fluidEnergySolver/fluidEnergySolver.H"
#include "solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.H"
#include "solverObjects/MULESVolumeFractionSolver/MULESVolumeFractionSolver.H"
#include "general/referenceFields/referenceFields.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchFields.H"
#include "fields/fvsPatchFields/constraint/empty/emptyFvsPatchField.H"
#include "primitives/functions/Function1/Table/TableBase.H"
#include "primitives/functions/Function1/Table/Table.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fv::flowSolver::correctInactiveGIBZoneFaces
(
    surfaceScalarField& fb
) const
{
    const volScalarField& p = this->p();
    const fvMesh& mesh = this->mesh();

    //- loop in the p field
    //  Check if there is an inactive GIB zone from BC
    //  If yes, grab cells and assign the at their faces fb = 0

    forAll(p.boundaryField(), pI)
    {
        const fvPatchScalarField& pfb = p.boundaryField()[pI];
        if (isA<fixedValueZoneFvPatchField<scalar>>(pfb))
        {
            const indirectPolyPatch& gibPolyPatch =
                refCast<const indirectPolyPatch>(pfb.patch().patch());

            const word czName =
                mesh.cellZones()[gibPolyPatch.zoneId()].name();

            const label& zoneId = mesh.cellZones().findZoneID(czName);
            const labelList& cz = mesh.cellZones()[zoneId];

            const cellList& cells = mesh.cells();
            forAll(cz, czI)
            {
                const label& gcI = cz[czI];
                forAll(cells[gcI], fI)
                {
                    const label gfI = cells[gcI][fI];
                    if (gfI< mesh.nInternalFaces())
                    {
                        fb[gfI] = 0;
                    }
                    else
                    {
                        label patchi = mesh.boundaryMesh().whichPatch(gfI);
                        label lfI = gfI - mesh.boundaryMesh()[patchi].start();
                        if
                        (
                            !isA<emptyFvsPatchField<scalar>>
                            (
                                fb.boundaryField()[patchi]
                            )
                        )
                        {
                            fb.boundaryFieldRef()[patchi][lfI] = 0;
                        }
                    }
                }
            }
            fb.boundaryFieldRef()[pI] = 0;
        }
    }
}


bool Foam::fv::flowSolver::foundThermoUpdateSolver()
{
    bool foundThermoUpdateSolver = false;
    const optionList& fvOptionsList = this->fvOptions();
    forAll(fvOptionsList, i)
    {
        if (fvOptionsList[i].type() == fluidEnergySolver::typeName)
        {
            foundThermoUpdateSolver = true;
        }
        if (fvOptionsList[i].type() == speciesConcentrationSolver::typeName)
        {
            foundThermoUpdateSolver = true;
        }
        if (fvOptionsList[i].type() == MULESVolumeFractionSolver::typeName)
        {
            foundThermoUpdateSolver = true;
        }
    }
    return foundThermoUpdateSolver;
}


void Foam::fv::flowSolver::makeTRefConstant()
{
    // Assume no changes in temperature
    referenceFields<scalar>& TRef =
        obr_.lookupObjectRef<referenceFields<scalar>>("TRef");
    // Value isn't important but good to avoid zero
    if (!TRef.isConst())
    {
        TRef.makeConst(300);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flowSolver::flowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    frameAcceleration_(nullptr),
    runTimeInfoDict_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flowSolver::initializeFrameAcceleration
(
    const dictionary& dict
)
{
    if (dict.found("frameAcceleration"))
    {
        frameAcceleration_.reset
        (
            Function1<vector>::New("frameAcceleration", dict).ptr()
        );
        createFrameAccelerationDict();
    }
}


void Foam::fv::flowSolver::initializeGravityHref()
{
    const Time& runTime = mesh().time();
    Info<< "Buoyancy active\n" << endl;
    Info<< "Reading g\n" << endl;
    g_.set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                obr(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "Reading hRef\n" << endl;
    hRef_.set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                obr(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimLength, 0)
        )
    );
    calculateghFields(true);
}


void Foam::fv::flowSolver::createFrameAccelerationDict()
{
    // Limiting function1 types to "constant" and "table"
    // The reason is to give to run-time post-processing
    // range (min/max). Which is easier from table and constant
    // functions.
    vector aVectorMax = Zero;
    vector aVectorMin = Zero;
    scalar aScalarMax = 0;
    scalar maxAPlusG = 0;

    // ensure that g_ is initialised
    dimensionedVector g("g", dimAcceleration, vector::zero);
    if (g_.valid())
    {
        g = g_();
    }

    const word Func1Type(frameAcceleration_().type());

    const Time& runTime = mesh().time();

    if (Func1Type == Function1Types::Table<vector>::typeName)
    {
        const Function1Types::Table<vector>& table =
            dynamic_cast<Function1Types::Table<vector>&>(frameAcceleration_());

        const vectorField y(table.y());

        const scalarField xComp(y.component(vector::X)());
        const scalarField yComp(y.component(vector::Y)());
        const scalarField zComp(y.component(vector::Z)());

        aVectorMax =
            vector
            (
                y[findMax(xComp)].x(),
                y[findMax(yComp)].y(),
                y[findMax(zComp)].z()
            );

        aVectorMin =
            vector
            (
                y[findMin(xComp)].x(),
                y[findMin(yComp)].y(),
                y[findMin(zComp)].z()
            );

        aScalarMax = mag(y[findMax(mag(y)())]);

        const scalarField magaPlusG(mag(y - g.value()));

        maxAPlusG = magaPlusG[findMax(magaPlusG)];
    }
    else if (Func1Type != Function1Types::Constant<vector>::typeName)
    {
        const vector constValue
        (
            frameAcceleration_->value(runTime.timeOutputValue())
        );
        aVectorMax = constValue;
        aVectorMin = constValue;
        aScalarMax = mag(constValue);
        maxAPlusG = mag(constValue - g.value());
    }
    else
    {
        FatalErrorInFunction
            << "frameAcceleration doesn't support: " << nl
            << Func1Type
            << " Function1 type." << nl
            << abort(FatalError);
    }

    // scalar sub-dict
    dictionary dictScalar("scalar");
    dictScalar.add("maxFrameAccelerationMagnitude", aScalarMax);
    dictScalar.add
    (
        "maxFrameAccelerationMagnitudeMinusGravity",
        maxAPlusG
    );

    // vector sub-dict
    dictionary dictVector("vector");
    dictVector.add
    (
        "frameAcceleration",
        frameAcceleration_->value(runTime.timeOutputValue())
    );
    dictVector.add("g", g.value());
    dictVector.add("maxFrameAccelerationComponents", aVectorMax);
    dictVector.add("minFrameAccelerationComponents", aVectorMin);

    // flowSolver sub-dict
    dictionary dictFlowSolver("flowSolver");
    dictFlowSolver.add("vector", dictVector);
    dictFlowSolver.add("scalar", dictScalar);

    // Results sub-dict
    dictionary dictAcceleration("results");
    dictAcceleration.add("flowSolver", dictFlowSolver);

    // One more dummy level to store correctly into the IOdictionary
    dictionary dictResults("results");
    dictResults.add("results", dictAcceleration);
    if (runTime.found("runTimeInfo"))
    {
        runTimeInfoDict_.set
        (
            runTime.lookupObjectRefPtr<IOdictionary>("runTimeInfo")
        );
        if (!runTimeInfoDict_().found("results"))
        {
            runTimeInfoDict_().set("results", dictionary("results"));
        }
        runTimeInfoDict_().subDict("results").set
        (
            "flowSolver",
            dictFlowSolver
        );
    }
    else
    {
        // Should initialise only when it is not available.
        runTimeInfoDict_.set
        (
            new IOdictionary
            (
                IOobject
                (
                    "runTimeInfo",
                    runTime.timeName(),
                    "uniform",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                dictResults
            )
        );
    }
}


void Foam::fv::flowSolver::calculateghFields(bool buoyant)
{
    if (!buoyant)
    {
        return;
    }

    const dimensionedScalar ghRef
    (
        mag(g_->value()) > SMALL
      ? g_() & (cmptMag(g_->value())/mag(g_->value()))*hRef_()
      : dimensionedScalar("ghRef", g_->dimensions()*dimLength, 0)
    );
    if (!gh_.valid())
    {
        Info<< "Creating field gh\n" << endl;
        gh_.set(new volScalarField("gh", (g_() & mesh().C()) - ghRef));
        ghf_.set(new surfaceScalarField("ghf", (g_() & mesh().Cf()) - ghRef));
    }
    else
    {
        gh_() = (g_() & mesh().C()) - ghRef;
        ghf_() = (g_() & mesh().Cf()) - ghRef;
    }

    // Correct boundary conditions to sync coupled boundaries as mesh.C() is a
    // sliced field and is not consistent in that regard
    gh_().correctBoundaryConditions();

    // Non-orthogonal correction for boundary gradient
    forAll(gh_->boundaryField(), patchi)
    {
        const vectorField pd(mesh().boundary()[patchi].delta());
        const labelUList& fc = mesh().boundary()[patchi].faceCells();
        if
        (
            !isA<emptyFvPatchField<scalar>>(gh_->boundaryField()[patchi])
            && !gh_->boundaryField()[patchi].coupled()
        )
        {
            fvPatchScalarField& pgh = gh_->boundaryFieldRef()[patchi];
            forAll(pgh, bfi)
            {
                pgh[bfi] =
                    (g_->value() & (mesh().C()[fc[bfi]] + pd[bfi]))
                    - ghRef.value();
            }
            ghf_->boundaryFieldRef()[patchi].forceAssign
            (
                gh_->boundaryField()[patchi]
            );
        }
    }
}


void Foam::fv::flowSolver::calculateFrameAccelerationContribution()
{
    if (frameAcceleration_.valid())
    {
        const dimensionedVector fAcc
        (
            "fAcc",
            dimAcceleration,
            frameAcceleration_->value(mesh_.time().value())
        );

        if (!gh_.valid())
        {
            Info<< "Creating field gh\n" << endl;
            gh_.set(new volScalarField("gh", (-fAcc & mesh().C())));
            ghf_.set(new surfaceScalarField("ghf", (-fAcc & mesh().Cf())));
        }
        else
        {
            gh_() += (-fAcc & mesh().C());
            ghf_() += (-fAcc & mesh().Cf());
        }

        // Updating RTPP dict
        runTimeInfoDict_().subDict("results").subDict("flowSolver").subDict("vector").set
        (
            "frameAcceleration",
            fAcc.value()
        );
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fv::flowSolver::faceBuoyancyForce
(
    const word& pDiffName,
    bool includeSnGradP
) const
{
    tmp<volScalarField> tbRho = buoyantRho();
    const volScalarField& bRho = tbRho();

    // Add and subtract snGrad p so that (p-rho*gh) is lumped into a single
    // snGrad, to minimise clipping by limited scheme
    // Ensure that we use the same snGrad scheme as the pressure laplacian,
    // for consistency
    tmp<fv::laplacianScheme<scalar, scalar>> pLaplacianScheme
    (
        fv::laplacianScheme<scalar, scalar>::New
        (
            mesh_,
            p().db(),
            mesh_.schemes().laplacianScheme
            (
                "laplacian(" + pDiffName + "," + p().name() + ')'
            ),
            "grad(" + p().name() + ')'
        )
    );
    tmp<surfaceScalarField> fbSource =
        - ghf_()*pLaplacianScheme().snGradientScheme().snGrad(bRho)
        - pLaplacianScheme().snGradientScheme().snGrad(p() - bRho*gh_());
    // Note: if includeSnGradP==true, this leaves it unbalanced, thereby
    // including it
    if (!includeSnGradP)
    {
        fbSource.ref() += pLaplacianScheme().snGradientScheme().snGrad(p());
    }
    correctInactiveGIBZoneFaces(fbSource.ref());

    return fbSource;
}


void Foam::fv::flowSolver::setOrComputeRhof()
{
    if (!thermo().isochoric() || !rhof_.valid() || !isStatic())
    {
        rhof_.clear();
        rhof_.reset(new surfaceScalarField("rhof", rhoInterpolation()));
    }
}


const Foam::tmp<Foam::surfaceScalarField>
Foam::fv::flowSolver::rhoInterpolation(const volScalarField& rho) const
{
    if (solnControl().transonic())
    {
        return
            fvc::interpolate
            (
                rho,
                phiv(),
                "interpolate(" + phiv().name() + "," + rho.name() + ")"
            );
    }
    else
    {
        return fvc::interpolate(rho);
    }
}


// ************************************************************************* //
