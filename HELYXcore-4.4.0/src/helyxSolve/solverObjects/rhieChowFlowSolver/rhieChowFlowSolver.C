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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rhieChowFlowSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/CourantNumber/CourantNumber.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "cfdTools/general/adjustSplitBCs/adjustSplitBCs.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"
#include "cfdTools/general/limitVolumes/limitVolumes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(rhieChowFlowSolver, 0);
}
}

makeFvSolverOption(rhieChowFlowSolver);

const word Foam::fv::rhieChowFlowSolver::rhorAUfName_ = "rhorAUf";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::rhieChowFlowSolver::rhieChowFlowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    flowSolver(name, obr, dict),
    thermoPtr_(nullptr),
    ddtPhiCorr_(true),
    printContErr_(true),
    solnControlPtr_(nullptr),
    useGradP_(false),
    transient_(false),
    buoyant_(false),
    distinctBuoyancy_(false),
    needsInitPhi_(false),
    needsCorrectPhi_(false)
{
    const Time& runTime = mesh().time();

    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    transient_ = !isA<simpleControl>(*solnControlPtr_);

    if (!obr_.foundObject<volScalarField>("p"))
    {
        Info<< "Creating initial p field\n" << endl;
        p_ = new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );
        obr_.store(p_);
    }
    else
    {
        Info<< "Reading field p\n" << endl;
        p_ = obr.lookupObjectRefPtr<volScalarField>("p");
    }
    volScalarField& p(*p_);

    Info<< "Reading field U\n" << endl;
    U_.set
    (
        new volVectorField
        (
            IOobject
            (
                IOobject::groupName("U", phaseName_),
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<rhoThermo>(multiphaseThermo::lookupOrCreate(obr, phaseName_));

    // Initialize info from material models
    buoyant_ = thermo().buoyant();
    distinctBuoyancy_ = thermo().distinctBuoyancy();

    Info<< "Creating field rho\n" << endl;
    rho_.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("rho", phaseName_),
                runTime.timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo().rho()
        )
    );

    if (buoyant_)
    {
        Info<< "Creating buoyant rho\n" << endl;
        bRho_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("buoyantRho", phaseName_),
                    runTime.timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_()
            )
        );
        if (distinctBuoyancy_)
        {
            bRho_() = thermo().buoyantRho();
        }

        initializeGravityHref();
    }
    initializeFrameAcceleration(dict);
    calculateFrameAccelerationContribution();

    updateSettingsFromDict();

    pressureControl_.set
    (
        new pressureControl
        (
            p,
            rho_,
            solnControlPtr_->dict(),
            thermo().pRefValue(),
            p.needReference() && (!transient_ || thermo().incompressible())
        )
    );

    if (!thermo().incompressible())
    {
        compressibility_ =
            new dimensionedScalar(fvc::domainIntegrate(thermo().psi()));
    }
    if (p.needReference() && thermo().incompressible())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureControl_->refValue()
          - getRefCellValue(p, pressureControl_->refCell())
        );
    }

    if (verbose_)
    {
        Info<< "p min/max: "
            << min(p).value() << " "
            << max(p).value() << endl;
        Info<< "rho min/max: "
            << min(rho_()).value() << " "
            << max(rho_()).value() << endl;

        if (distinctBuoyancy_)
        {
            Info<< "Buoyant rho min/max: "
                << min(bRho_()).value() << " "
                << max(bRho_()).value() << endl;
        }
    }

    mesh().schemes().setFluxRequired(p.name());

    initialMass_.set(new dimensionedScalar(fvc::domainIntegrate(rho_())));

    rhoMin_.set
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMin",
                solnControlPtr_->dict(),
                dimDensity,
                0
            )
        )
    );

    rhoMax_.set
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMax",
                solnControlPtr_->dict(),
                dimDensity,
                GREAT
            )
        )
    );

    // Read or create face flux
    IOobject phiIOObj
    (
        IOobject::groupName("phi", phaseName_),
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phiIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating face flux field phi" << nl << endl;
        needsInitPhi_ = true;
    }
    else
    {
        Info<< "Reading face flux field phi" << nl << endl;
        needsInitPhi_ = false;
    }
    phi_.set
    (
        new surfaceScalarField
        (
            phiIOObj,
            linearInterpolate(rho_()*U_()) & mesh().Sf()
        )
    );

    if (solnControlPtr_->transonic())
    {
        phiv_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phiv", phaseName_),
                    runTime.timeName(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(U_()) & mesh().Sf()
            )
        );
        // Make flux consistent with density interpolation
        phi_() = rhoInterpolation()*phiv_();
    }

    if (!isStatic())
    {
        fvc::makeRelative(phi_(), rho_(), U_());
    }

    cumulativeContErr_ = 0;

    Info<< "Creating turbulence model\n" << endl;
    turbulence_=
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            phi_,
            thermo()
        ).ptr();
    turbulence_->store();
    limitVolumes(this->dict(), mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rhieChowFlowSolver::reactingStoreDpdt()
{
    if (isReactingFoam() && thermoPtr_->dpdt())
    {
        volScalarField::New
        (
            "dpdt",
            obr_,
            mesh_,
            dimensionedScalar(dimPressure/dimTime, 0)
        ).ptr()->store();
    }
}


void Foam::fv::rhieChowFlowSolver::reactingCorrect()
{
    if (isReactingFoam())
    {
        if (obr_.foundObject<volScalarField>("K"))
        {
            obr_.lookupObjectRef<volScalarField>("K") = 0.5*magSqr(U_());
        }
        if (pressureControl_->limit(*p_))
        {
            (*p_).correctBoundaryConditions();
        }
        rho_() = thermo().rho();
        if (thermoPtr_->dpdt())
        {
            obr_.lookupObjectRef<volScalarField>("dpdt") = fvc::ddt(*p_);
        }
    }
}


void Foam::fv::rhieChowFlowSolver::updateSettingsFromDict()
{
    ddtPhiCorr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "ddtPhiCorr",
        true
    );

    printContErr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "printContinuityErrors",
        true
    );

    useGradP_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "useGradP",
        !buoyant_
    );
}


void Foam::fv::rhieChowFlowSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    word UPredName = U_->name() + "Predictor";
    solveNames = {"correctPhi", UPredName, p_->name(), U_->name()};

    // Dependencies
    optionalDependencies.insert(UPredName, {"fvMesh"});
    requiredDependencies.insert(p_->name(), {UPredName});
    requiredDependencies.insert(U_->name(), {p_->name()});

    optionalDependencies.insert("correctPhi", {"fvMesh"});
    requiredDependencies.insert(UPredName, {"correctPhi"});

    optionalDependencies.insert("turbulence", solveNames);

    // Correctors

    // Add all solves to outer corrector
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);

    // PISO (inner) correctors
    correctorMembers.insert("PISOCorrector", {p_->name(), U_->name()});

    // Non-orthogonal correctors
    correctorMembers.insert
    (
        "nonOrthogonalCorrector:" + p_->name(), {p_->name()}
    );

    // Switching solution order
    if (isReactingFoam())
    {
        optionalDependencies.insert(p_->name(), {thermoPtr_->heT().name()});
    }

    // optional initialisation step: correctPhi after mesh update
    solveNames.append("solverMeshInit");
    optionalDependencies.insert("solverMeshInit", {"fvMeshInit"});
    correctorMembers.insert
    (
        solverObject::initialisationLoopName, {"solverMeshInit"}
    );
}


bool Foam::fv::rhieChowFlowSolver::initialise()
{
    if (transient_ && (!isStatic() || fvOptions().hasMRF()))
    {
        Info<< "Reading/calculating field rhoUf\n" << endl;
        rhoUf_.set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName("rhoUf", phaseName_),
                    mesh().time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho_()*U_())
            )
        );
    }
    volScalarField& p = *p_;

    if (needsInitPhi_ && !transient_)
    {
        // Set fixed-flux pressure gradient to its existing value simple to
        // avoid a warning in matrix construction for fvOptions below.
        // We do not want to set it to zero in case field initialisation has
        // already assigned it a correct value.
        volScalarField::Boundary& pBf = p.boundaryFieldRef();
        forAll(pBf, patchi)
        {
            if (isA<fixedFluxPressureFvPatchScalarField>(pBf[patchi]))
            {
                fixedFluxPressureFvPatchScalarField& ffp =
                    refCast<fixedFluxPressureFvPatchScalarField>
                    (
                        pBf[patchi]
                    );
                ffp.updateSnGrad(scalarField(ffp.gradient()));
            }
        }

        tmp<volScalarField> fvOptSrc
        (
            fvOptions()(thermo().psi(), p, rho_().name()) & p
        );

        CorrectPhi
        (
            U_(),
            phi_(),
            p,
            rho_(),
            thermo().psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            fvOptSrc(),
            *solnControlPtr_,
            p.name(),
            (p.needReference() && (!transient_ || thermo().incompressible())),
            pressureControl_->refCell()
        );
        // updateCoeffs was called in matrix construction but evaluate was
        // not called - resetUpdate reverts to updatable state
        forAll(p.boundaryField(), patchi)
        {
            p.boundaryFieldRef()[patchi].resetUpdate();
        }

        Info<< endl;
    }

    if (!exists(phi_().objectPath()))
    {
        tmp<surfaceScalarField> rhof(fvc::interpolate(rho_()));
        fvOptions().makeRelative(rhof, phi_());
    }

    if (transient_)
    {
        // Calculate initial Courant number
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());
    }
    reactingStoreDpdt();

    return true;
}


scalar Foam::fv::rhieChowFlowSolver::getMaxTimeStep()
{
    // Try region-specific maxCo, otherwise revert to global
    const dictionary& maxCoDict
    (
        solnControlPtr_->dict().found("maxCo")
      ? solnControlPtr_->dict()
      : mesh_.time().controlDict()
    );
    maxCoDataPtr_ = Function1<scalar>::New("maxCo", maxCoDict);

    scalar maxMaxCo = maxCoDataPtr_->value(mesh_.time().timeIndex());
    return maxMaxCo/stabilise(maxCo_, SMALL)*mesh_.time().deltaTValue();
}


bool Foam::fv::rhieChowFlowSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if (correctorName == "PISOCorrector")
    {
        if (transient_)
        {
            label nCorr = refCast<pimpleControl>(*solnControlPtr_).nCorrPISO();
            return (corrector >= nCorr-1);
        }
        else
        {
            return true;
        }
    }
    else if (correctorName == "nonOrthogonalCorrector:" + p_->name())
    {
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return true;
    }
}


void Foam::fv::rhieChowFlowSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    // Update every time step
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        updateSettingsFromDict();
    }

    if (correctorName == solverObject::outerCorrectorName)
    {
        if (!isStatic())
        {
            // Store divrhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
            divRhoU_.clear();
            divRhoU_ =
                tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "divrhoU",
                        fvc::div
                        (
                            fvOptions().absolute
                            (
                                fvc::absolute(phi_(), rho_(), U_()),
                                rho_()
                            )
                        )
                    )
                );
        }
    }
    else if (correctorName == "PISOCorrector")
    {
        corr_ = corrector;
    }
    else if (correctorName == "nonOrthogonalCorrector:" + p_->name())
    {
        nonOrthCorr_ = corrector;
    }
}


tmp<fvVectorMatrix>
Foam::fv::rhieChowFlowSolver::assembleVectorMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if (fieldName == U_->name() + "Predictor")
    {
        if (correctorNumber_[solverObject::outerCorrectorName] == 0)
        {
            // Explicit predictor for new value of rho
            // Only use if compressible and running in PISO mode
            if
            (
                transient_
             && !thermo().isochoric()
             && solnControlPtr_->nOuterCorr() <= 1
            )
            {
                solveContinuity();
            }
            rho_().storePrevIter();
        }
        if (solnControlPtr_->momentumPredictor())
        {
            return UPredictorEqn();
        }
        else
        {
            return tmp<fvVectorMatrix>();
        }
    }
    else
    {
        return tmp<fvVectorMatrix>();
    }
}


tmp<fvVectorMatrix> Foam::fv::rhieChowFlowSolver::UPredictorEqn()
{
    assembleUEqnLHS();

    tmp<fvVectorMatrix> tUEqn(new fvVectorMatrix(UEqnLHS_()));
    if (useGradP_)
    {
        tUEqn.ref() += fvc::grad(*p_);
        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce(rhorAUfName_);
            tUEqn.ref() -=
                fvc::reconstruct(faceBuoyancyForce_()*mesh_.magSf());
        }
    }
    else
    {
        tmp<surfaceScalarField> tsnGrad;
        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce(rhorAUfName_, true);
            tsnGrad = -faceBuoyancyForce_();
        }
        else
        {
            tsnGrad = fvc::snGrad(*p_);
        }
        tUEqn.ref() += fvc::reconstruct(tsnGrad*mesh_.magSf());
    }
    return tUEqn;
}


void Foam::fv::rhieChowFlowSolver::assembleUEqnLHS()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();

    fv::options& fvOptions = this->fvOptions();

    if (!transient_ && solnControlPtr_->modifiedMomentumInterp())
    {
        U.storePrevIter();
    }

    // Assemble the Momentum equation
    UEqnLHS_ =
        new fvVectorMatrix
        (
            fvm::div(phi, U)
          + fvOptions.MRFDDt(rho, U)
          + turbulence_->divDevRhoReff(U)
        ==
            fvOptions(rho, U)
        );

    if (transient_)
    {
        if (!thermo().isochoric() || !isStatic())
        {
            UEqnLHS_.ref() += fvm::ddt(rho, U);
        }
        else
        {
            const word schemeName("ddt(" + rho.name() + "," + U.name() + ")");
            UEqnLHS_.ref() += rho*fvm::ddt(U, schemeName);
        }
    }

    UEqnLHS_->relax();

    fvOptions.constrain(UEqnLHS_.ref());
}


tmp<fvScalarMatrix>
Foam::fv::rhieChowFlowSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if (fieldName == p_->name())
    {
        assemblepEqn();
        finalSolve =
            (
                transient_
             && finalIter_[solverObject::outerCorrectorName]
             && finalIter_["PISOCorrector"]
             && finalIter_["nonOrthogonalCorrector:" + p_->name()]
            );
        // Return copy of the tmp so it doesn't get deleted from outside
        return tmp<fvScalarMatrix>(pEqn_);
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }
}


void Foam::fv::rhieChowFlowSolver::assemblepEqn()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();
    volScalarField& p = *p_;
    const volScalarField& psi = thermo().psi();
    fv::options& fvOptions = this->fvOptions();

    if (frameAcceleration_.valid())
    {
        calculateghFields(buoyant_);
    }
    calculateFrameAccelerationContribution();

    if
    (
        !solnControlPtr_->momentumPredictor()
     && corr_ == 0
     && nonOrthCorr_ == 0
    )
    {
        assembleUEqnLHS();
    }

    if (nonOrthCorr_ == 0)
    {
        if (!thermo().isochoric())
        {
            rho = thermo().rho();

            // Relax density here in response to temperature change. Relaxed again
            // in correct() to respond to pressure change
            // Store for consistency in case no relaxation in use
            rho.relax();
        }
        if (buoyant_)
        {
            if (distinctBuoyancy_)
            {
                bRho_() = thermo().buoyantRho();
            }
            else
            {
                bRho_() = rho;
            }
        }

        if (!thermo().incompressible())
        {
            // Thermodynamic density needs to be updated by psi*d(p) after the
            // pressure solution
            p0_.clear();
            p0_ = volScalarField::New(typeName + "-p0", p);
            p0_.ref().checkIn();
        }

        // Construct unregistred and decide about registration later
        tmp<volScalarField> trAU
        (
            volScalarField::New("(1|A(U))", 1/UEqnLHS_().A())
        );

        HbyA_ = constrainHbyA(trAU()*UEqnLHS_().H(), U, p);
        volVectorField& HbyA = HbyA_.ref();
        if (solnControl().consistent())
        {
            tmp<volScalarField> A(1/trAU());
            tmp<volScalarField> AmH1(-UEqnLHS_().H1() + A());
            // Stabilise the case where H1 ~= -A
            rAtU_ = volScalarField::New("(1|At(U))", 1/max(AmH1, 0.01*A));
            trAU.ref().checkIn();
        }
        else
        {
            // Transfer the tmp
            rAtU_ = trAU;
            rAtU_.ref().checkIn();
        }
        const volScalarField& rAtU = rAtU_();

        tmp<surfaceScalarField> phivHbyA
        (
            surfaceScalarField::New("phivHbyA", fvc::flux(HbyA))
        );
        tmp<surfaceScalarField> rAtUf
        (
            surfaceScalarField::New
            (
                "rAUf",
                fvc::interpolate(rAtU, "interpolate((1|A(U)))")
            )
        );

        setOrComputeRhof();

        if (solnControlPtr_->transonic())
        {
            rhorAtUf_ = surfaceScalarField::New(rhorAUfName_, rhof_()*rAtUf());
            phiHbyA_ = surfaceScalarField::New("phiHbyA", rhof_()*phivHbyA());
        }
        else
        {
            // Reuse tmps for efficiency
            rhorAtUf_ = surfaceScalarField::New(rhorAUfName_, rhof_()*rAtUf);
            phiHbyA_ = surfaceScalarField::New("phiHbyA", rhof_()*phivHbyA);
        }
        const surfaceScalarField& rhorAtUf = rhorAtUf_();
        surfaceScalarField& phiHbyA = phiHbyA_.ref();

        tmp<surfaceScalarField> rhorAUf;
        if (solnControlPtr_->consistent())
        {
            tmp<surfaceScalarField> rAUf =
                fvc::interpolate(trAU(), "interpolate((1|A(U)))");
            rhorAUf = rhof_()*rAUf;
        }
        else
        {
            // Point tmp to existing var; increments ref count
            rhorAUf = tmp<surfaceScalarField>(rhorAtUf_);
        }

        if (transient_ && ddtPhiCorr_)
        {
            phiHbyA += rhorAUf()*fvc::ddtCorr(rho, U, phi, rhoUf_);
        }
        fvOptions.makeRelative(rhof_(), phiHbyA);

        // If phi was (re)initialised with correctPhi, this causes a large
        // discrepancy here between phi and U, so we skip it
        if
        (
            solnControlPtr_->modifiedMomentumInterp()
         && !transient_
         && !needsCorrectPhi_ && !needsInitPhi_
        )
        {
            tmp<surfaceScalarField> rho0f =
                fvc::interpolate(rho.prevIter(), "interpolate(rho)");
            fvOptions.makeAbsolute(rho0f(), phi);
            phiHbyA +=
                (1 - mesh_.solution().equationRelaxationFactor(U.name()))
               *rho0f()
               *(
                    fvc::absolute(phi/rho0f(), U)
                  - (fvc::interpolate(U.prevIter()) & mesh_.Sf())
                );
        }
        needsInitPhi_ = false;
        needsCorrectPhi_ = false;

        if (transient_)
        {
            closedVolume_ = p.needReference();
        }
        else
        {
            closedVolume_ = adjustPhi(phiHbyA, U, p);
        }

        if (solnControl().consistent())
        {
            tmp<surfaceScalarField> tsnGrad;
            if (useGradP_)
            {
                tsnGrad = fvc::snGrad(p)*mesh().magSf();
                tmp<volVectorField> tgrad(fvc::grad(p));
                if (buoyant_)
                {
                    tsnGrad.ref() -= faceBuoyancyForce_()*mesh_.magSf();
                    tgrad.ref() -=
                        fvc::reconstruct(faceBuoyancyForce_*mesh_.magSf());
                }
                phiHbyA += (rhorAtUf-rhorAUf())*tsnGrad();
                HbyA -= (trAU() - rAtU)*tgrad;
            }
            else
            {
                if (buoyant_)
                {
                    tsnGrad = -faceBuoyancyForce_*mesh().magSf();
                }
                else
                {
                    tsnGrad = fvc::snGrad(p)*mesh().magSf();
                }
                phiHbyA += (rhorAtUf-rhorAUf())*tsnGrad();
                HbyA -= (trAU() - rAtU)*fvc::reconstruct(tsnGrad);
            }
        }

        if (buoyant_)
        {
            faceBuoyancyForce_ = this->faceBuoyancyForce(rhorAUfName_);
            phiHbyA += rhorAtUf_()*faceBuoyancyForce_()*mesh_.magSf();
        }

        // Update the fixedFluxPressure BCs to ensure flux consistency
        constrainPressure(p, rho, U, phiHbyA, rhorAtUf_(), fvOptions);

        if (!thermo().incompressible())
        {
            compressibility_() = fvc::domainIntegrate(psi);
        }

        // For correct behaviour this term has to be after constrain pressure
        // It compensates the fvc::ddt(rho) that produce additional flux (mesh)
        // Even for constant rho
        if (!isStatic())
        {
            fvc::makeRelative(phiHbyA, rho, U);
        }

        pDDtEqn_ =
            new fvScalarMatrix
            (
              - fvOptions(psi, p, rho.name())
              + fvc::div(phiHbyA)
            );
        if (!thermo().isochoric() || !isStatic())
        {
            pDDtEqn_.ref() += fvc::ddt(rho);
        }
        if (!thermo().incompressible())
        {
            pDDtEqn_.ref() += psi*correction(fvm::ddt(p));
        }

        if (solnControlPtr_->transonic())
        {
            phiv_() =
                phivHbyA - rAtUf*fvc::snGrad(p)*mesh_.magSf();
            gaussConvectionScheme<scalar> convScheme
            (
                mesh_,
                phiv_(),
                tmp<surfaceInterpolationScheme<scalar>>
                (
                    new upwind<scalar>(mesh_, phiv_())
                )
            );
            pDDtEqn_.ref() += correction(convScheme.fvmDiv(phiv_(), psi, p));
        }
    }

    // Solve pressure
    pEqn_ =
        new fvScalarMatrix(pDDtEqn_() - fvm::laplacian(rhorAtUf_(), p));

    if (solnControlPtr_->transonic())
    {
        pEqn_->relax();
    }

    // Do not set reference in transient, compressible case as there is a dp/dt
    // term
    if (p.needReference() && (!transient_ || thermo().incompressible()))
    {
        // In the steady, compressible case do not reset to ref. value as the
        // pressure is later modified to respect global mass conservation
        pEqn_->setReference
        (
            pressureControl_->refCell(),
            thermo().incompressible()
          ? pressureControl_->refValue()
          : getRefCellValue(p, pressureControl_->refCell())
        );
    }
}


void Foam::fv::rhieChowFlowSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();
    const word pCorrName("nonOrthogonalCorrector:" + p_->name());
    if (solveName == U_->name() + "Predictor")
    {
        fvOptions.correct(U_());

        if (verbose_)
        {
            Info<< "U predictor: mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }
    }
    else if (solveName == p_->name())
    {
        if (verbose_)
        {
            Info<< "p min/max: "
                << min(*p_).value() << " "
                << max(*p_).value() << endl;
        }

        if (finalIter_[pCorrName])
        {
            // Final non-orthogonal corrector

            volScalarField& p = *p_;
            volScalarField& rho = rho_();

            // Calculate the conservative flux
            phi_() = phiHbyA_ + pEqn_->flux();

            adjustSplitBCs(phi_(), U_(), p);
            adjustSplitBCs(phi_(), U_(), U_());
            p.relax();
            pressureControl_->limit(p);

            if (printContErr_)
            {
                if (!thermo().isochoric())
                {
                    if (transient_)
                    {
                        // Here this is only used for printing continuity errors
                        // as rho is overwritten below
                        solveContinuity();
                    }
                    printCompContinuityErrors(regionName);
                }
                else
                {
                    printIncoContinuityErrors(regionName);
                }
            }

            if (closedVolume_)
            {
                if (thermo().incompressible())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pressureControl_->refValue()
                      - getRefCellValue(p, pressureControl_->refCell())
                    );
                }
                else
                {
                    // Ensure that density correction conserves
                    // initial mass
                    p +=
                        (
                            initialMass_()
                          - fvc::domainIntegrate
                            (
                                thermo().rho() + thermo().psi()*(p-p0_())
                            )
                        )/compressibility_();
                }
            }

            if (!thermo().incompressible())
            {
                thermoPtr_->correctRho(thermo().psi()*(p - p0_()));
            }

            if (!thermo().isochoric())
            {
                rho = thermo().rho();
                if (!transient_)
                {
                    rho = max(rho, rhoMin_());
                    rho = min(rho, rhoMax_());
                }
                rho.relax();
            }
            else
            {
                // Keep time index up to date
                rho.timeIndex() = mesh().time().timeIndex();
            }
        }
    }
    if
    (
        (solveName == U_->name() && !isReactingFoam())
     || (
           isReactingFoam() && solveName == p_->name() && finalIter_[pCorrName]
        )
    )
    {
        // U corrector

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        volVectorField& U = U_();
        volScalarField& p = *p_;
        volScalarField& rho = rho_();
        surfaceScalarField& phi = phi_();
        if (useGradP_)
        {
            U = HbyA_ - rAtU_()*fvc::grad(p);
            if (buoyant_)
            {
                U += rAtU_()*fvc::reconstruct(faceBuoyancyForce_*mesh_.magSf());
            }
        }
        else
        {
            //tmp<surfaceScalarField> pEqnFlux(pEqn_->flux()/rhorAtUf_());
            // At present we can't use flux here because pressure relaxation
            // may affect pressure boundaries, yet this is frozen into the
            // boundaryCoeffs so will not reflect in flux
            tmp<surfaceScalarField> pEqnFlux = -fvc::snGrad(p);
            if (buoyant_)
            {
                pEqnFlux.ref() += faceBuoyancyForce_();
            }
            pEqnFlux.ref() *= mesh_.magSf();
            U = HbyA_ + rAtU_()*fvc::reconstruct(pEqnFlux);
        }

        pEqn_.clear();

        U.correctBoundaryConditions();
        fvOptions.correct(U);

        if (verbose_)
        {
            Info<< "mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }

        if (rhoUf_.valid())
        {
            rhoUf_() = fvc::interpolate(rho*U);
            surfaceVectorField n(mesh_.Sf()/mesh_.magSf());
            rhoUf_() +=
                n
                *(
                    fvOptions.absolute(fvc::absolute(phi, rho, U), rho)
                   /mesh_.magSf()
                  - (n & rhoUf_())
                );
        }

        if (finalIter_["PISOCorrector"] && transient_)
        {
            if (!thermo().isochoric())
            {
                // To undo the relaxation?
                rho = thermo().rho();
            }
        }
        reactingCorrect();
    }
    else if (solveName == "correctPhi")
    {
        if (needsCorrectPhi_)
        {
            correctPhi(true);
            calculateghFields(buoyant_);
            calculateFrameAccelerationContribution();
            // needsCorrectPhi_ is reset to false in
            // pressure equation construction
        }
    }
    else if
    (
        solveName == "solverMeshInit"
     && dict().lookupOrDefault<bool>("solverMeshInit", false)
    )
    {
        Info<< "Solve solverMeshInit" << endl;
        // Calculate absolute flux from the mapped surface velocity
        phi_() = mesh_.Sf() & fvc::interpolate(rho_()*U_());
        tmp<volScalarField> tdivU(fvc::div(phi_())*0.0);

        CorrectPhi
        (
            U_(),
            phi_(),
            *p_,
            rho_(),
            thermo().psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            tdivU(),
            *solnControlPtr_,
            word::null,
            (p_->needReference() && (!transient_ || thermo().incompressible())),
            pressureControl_->refCell()
        );

        // Make the fluxes relative to the mesh-motion
        fvc::makeRelative(phi_(), rho_(), U_());
    }

    // Recomputing material properties on the end of every outer iteration
    // when the fluid energy solver isn't available.
    // TODO: This is temporary solution. It should be handled more generally.
    // Maybe by putting material library update to solver object?
    if
    (
        finalIter_.found("PISOCorrector")
     && finalIter_["PISOCorrector"]
     && finalIter_["nonOrthogonalCorrector:p"]
     && solveName == U_->name()
    )
    {
        if (!foundThermoUpdateSolver() && !thermo().isConst())
        {
            makeTRefConstant();
            const_cast<fluidThermo&>(thermo()).correct();
        }
    }
}


void Foam::fv::rhieChowFlowSolver::solveContinuity()
{
    fv::options& fvOptions = this->fvOptions();
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_())
      + fvc::div(phi_())
      ==
        fvOptions(rho_())
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho_());
}


void Foam::fv::rhieChowFlowSolver::printIncoContinuityErrors
(
    const word& regionName
)
{
    volScalarField contErr(fvc::div(phi_()));

    const scalar sumLocalContErr =
        mesh_.time().deltaTValue()
       *mag(contErr)().weightedAverage(mesh_.V()).value();

    const scalar globalContErr =
        mesh_.time().deltaTValue()
       *contErr.weightedAverage(mesh_.V()).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::rhieChowFlowSolver::printCompContinuityErrors
(
    const word& regionName
)
{
    const volScalarField& rho = rho_();

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr, globalContErr;

    if (transient_)
    {
        sumLocalContErr =
            (
                fvc::domainIntegrate(mag(rho - thermo().rho()))/totalMass
            ).value();

        globalContErr =
            (
                fvc::domainIntegrate(rho - thermo().rho())/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;
    }
    else
    {
        dimensionedScalar totalMass = fvc::domainIntegrate(rho);

        fv::options& fvOptions = this->fvOptions();
        volScalarField rhoSource ( fvOptions(rho_()) & rho );
        sumLocalContErr =
            (
                fvc::domainIntegrate(mag(fvc::div(phi_()) - rhoSource))
               /totalMass
            ).value();

        globalContErr =
            (
                fvc::domainIntegrate(fvc::div(phi_())-rhoSource)/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;
    }

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::rhieChowFlowSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if
    (
        transient_
     && correctorName == solverObject::outerCorrectorName
     && finalIter
    )
    {
        Info<< endl;
        // Calculate Courant for next time step
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());
    }
}


void Foam::fv::rhieChowFlowSolver::correctPhi(bool defaultCorrectPhi)
{
    if
    (
        solnControlPtr_->dict().lookupOrDefault<Switch>
        (
            "correctPhi", defaultCorrectPhi
        )
    )
    {
        // Calculate absolute flux
        phi_() = mesh_.Sf() & fvc::interpolate(rho_()*U_());

        CorrectPhi
        (
            U_(),
            phi_(),
            *p_,
            rho_(),
            thermo().psi(),
            dimensionedScalar("rAUf", dimTime, 1),
            divRhoU_(),
            *solnControlPtr_,
            word::null,
            (p_->needReference() && (!transient_ || thermo().incompressible())),
            pressureControl_->refCell()
        );

        // Make the fluxes relative to the mesh-motion
        fvc::makeRelative(phi_(), rho_(), U_());
    }
}


bool Foam::fv::rhieChowFlowSolver::movePoints()
{
    if (mesh_.changing()) limitVolumes(this->dict(), mesh_);
    calculateghFields(buoyant_);
    calculateFrameAccelerationContribution();
    return true;
}


void Foam::fv::rhieChowFlowSolver::topoChange(const polyTopoChangeMap& map)
{
    if (mesh_.changing()) limitVolumes(this->dict(), mesh_);
    needsCorrectPhi_ = true;
}


bool Foam::fv::rhieChowFlowSolver::setRDeltaT()
{
    if (fv::localEulerDdt::enabled(mesh_))
    {
        const dictionary& solnDict = solnControlPtr_->dict();

        // The local time step will be registred under mesh object registry
        volScalarField& rDeltaT = lookupOrConstructRDeltaT();

        // Maximum flow Courant number
        const scalar maxCo = solnDict.lookup<scalar>("maxCo");

        tmp<volScalarField::Internal> rDeltaTT
        (
            volScalarField::Internal::New
            (
                typeName + "-rDeltaT",
                fvc::surfaceSum(mag(phi_()))()()/((2*maxCo)*mesh_.V()*rho_()())
            )
        );

        Info<< "    Flow = "
            << 1/max(gMax(rDeltaTT().field()), VSMALL) << ", "
            << 1/max(gMin(rDeltaTT().field()), VSMALL) << endl;

        rDeltaT.ref() = max(rDeltaT(), rDeltaTT());

        // Maximum time scale
        const scalar maxDeltaT =
            solnDict.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Limit max deltaT = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;

        // Update tho boundary values of the reciprocal time-step
        rDeltaT.correctBoundaryConditions();

        // Smoothing parameter (0-1) when smoothing iterations > 0
        const scalar rDeltaTSmoothingCoeff =
            solnDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1);

        // Spatially smooth the time scale field
        if (rDeltaTSmoothingCoeff < 1)
        {
            fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
        }

        // Limit rate of change of time scale
        // - reduce as much as required
        // - only increase at a fraction of old time scale

        // Damping coefficient (1-0)
        const scalar rDeltaTDampingCoeff =
            solnDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0);
        if
        (
            rDeltaTDampingCoeff < 1
         && mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1
        )
        {
            rDeltaT =
                max
                (
                    rDeltaT,
                    (scalar(1.0) - rDeltaTDampingCoeff)
                   *mesh_.lookupObject<volScalarField>("rDeltaT0")
                );
        }

        // Update tho boundary values of the reciprocal time-step
        rDeltaT.correctBoundaryConditions();

        Info<< "    Overall = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;

        return true;
    }

    return false;
}


// ************************************************************************* //
