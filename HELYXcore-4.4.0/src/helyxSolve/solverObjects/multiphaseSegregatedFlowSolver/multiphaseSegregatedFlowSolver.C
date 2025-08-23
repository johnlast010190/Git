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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseSegregatedFlowSolver.H"
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
#include "rhoThermo/rhoThermo.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "cfdTools/general/limitVolumes/limitVolumes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseSegregatedFlowSolver, 0);
}
}

makeFvSolverOption(multiphaseSegregatedFlowSolver);

const word Foam::fv::multiphaseSegregatedFlowSolver::rAUfName_ = "rAUf";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseSegregatedFlowSolver::multiphaseSegregatedFlowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    flowSolver(name, obr, dict),
    p_(nullptr),
    thermoPtr_(nullptr),
    closedVolume_(false),
    turbulence_(nullptr),
    printContErr_(true),
    cumulativeContErr_(0),
    maxCo_(0),
    solnControlPtr_(nullptr),
    buoyant_(false),
    needsCorrectPhi_(false),
    phiNeedsInit_(false),
    implicitUSolve_(false),
    implicitUCorrector_(false),
    momentumPredictor_(false),
    corr_(0),
    nonOrthCorr_(0)
{
    const Time& runTime = mesh().time();

    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    printContErr_ = solnControlPtr_->dict().lookupOrDefault<Switch>
    (
        "printContinuityErrors",
        true
    );

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
                "U",
                runTime.timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ = &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr));

    buoyant_ = thermoPtr_->buoyant();
    if (thermoPtr_->distinctBuoyancy())
    {
        FatalErrorInFunction
            << "This solver does not support distinct buoyancy." << nl << endl;
    }

    Info<< "Creating field rho\n" << endl;
    rho_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoPtr_->rho()
        )
    );

    if (buoyant_)
    {
        initializeGravityHref();
    }
    initializeFrameAcceleration(dict);
    calculateFrameAccelerationContribution();

    pressureControl_.set
    (
        new pressureControl
        (
            p,
            rho_,
            solnControlPtr_->dict(),
            thermoPtr_->pRefValue(),
            p.needReference() && thermoPtr_->incompressible()
        )
    );

    if (!thermoPtr_->incompressible())
    {
        compressibility_ =
            new dimensionedScalar(fvc::domainIntegrate(thermoPtr_->psi()));
    }
    if (p.needReference() && thermoPtr_->incompressible())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureControl_->refValue()
          - getRefCellValue(p, pressureControl_->refCell())
        );
    }

    mesh().schemes().setFluxRequired(p.name());

    initialMass_.set(new dimensionedScalar(fvc::domainIntegrate(rho_())));

    // Read or create volumetric face flux
    IOobject phivIOObj
    (
        "phiv",
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phivIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating volumetric face flux field phiv" << nl << endl;

        // We can't make phi relative until the mesh has been moved, so remember
        // that we need to do it later
        phiNeedsInit_ = true;
    }
    else
    {
        Info<< "Reading volumetric face flux field phiv" << nl << endl;
    }
    phiv_.set
    (
        new surfaceScalarField
        (
            phivIOObj,
            linearInterpolate(U_()) & mesh().Sf()
        )
    );

    // Read or create face flux
    IOobject phiIOObj
    (
        "phi",
        runTime.timeName(),
        obr,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!phiIOObj.typeHeaderOk<surfaceScalarField>())
    {
        Info<< "Calculating face flux field phi" << nl << endl;
    }
    else
    {
        Info<< "Reading face flux field phi" << nl << endl;
    }
    phi_.set
    (
        new surfaceScalarField
        (
            phiIOObj,
            linearInterpolate(rho_())*phiv_()
        )
    );

    cumulativeContErr_ = 0;

    Info<< "Creating turbulence model\n" << endl;
    turbulence_ =
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            phi_,
            *thermoPtr_
        ).ptr();
    turbulence_->store();
    limitVolumes(this->dict(), mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseSegregatedFlowSolver::getSolveGraph
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
    optionalDependencies.insert(UPredName, {"alphas"});
    requiredDependencies.insert(p_->name(), {UPredName});
    requiredDependencies.insert(U_->name(), {p_->name()});

    optionalDependencies.insert("correctPhi", {"fvMesh"});
    requiredDependencies.insert("alphas", {"correctPhi"});

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
}


bool Foam::fv::multiphaseSegregatedFlowSolver::initialise()
{
    if (!isStatic() || fvOptions().hasMRF())
    {
        Info<< "Reading/calculating field Uf\n" << endl;
        Uf_.set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    mesh().time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U_())
            )
        );
    }

    if (!fv::localEulerDdt::enabled(mesh_))
    {
        // Calculate initial Courant number
        maxCo_ =
            CourantNumber(mesh_, mesh_.time().deltaTValue(), rho_(), phi_());
    }

    return true;
}


scalar Foam::fv::multiphaseSegregatedFlowSolver::getMaxTimeStep()
{
    if (!fv::localEulerDdt::enabled(mesh_))
    {
        scalar maxMaxCo = 0;

        // Try region-specific maxCo, otherwise revert to global
        const dictionary& solnDict
        (
            solnControlPtr_->dict().found("maxCo")
        ? solnControlPtr_->dict()
        : mesh_.time().controlDict()
        );

        maxCoDataPtr_ = Function1<scalar>::New("maxCo", solnDict);

        maxMaxCo = maxCoDataPtr_->value(mesh_.time().timeIndex());

        return maxMaxCo/stabilise(maxCo_, SMALL)*mesh_.time().deltaTValue();
    }
    else
    {
        return GREAT;
    }
}


bool Foam::fv::multiphaseSegregatedFlowSolver::isFinalCorrector
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
        label nCorr = refCast<pimpleControl>(*solnControlPtr_).nCorrPISO();
        return (corrector >= nCorr-1);
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return false;
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        implicitUCorrector_ =
            solnControlPtr_->dict().lookupOrDefault
            (
                "implicitUCorrector", false
            );
        implicitUSolve_ =
            solnControlPtr_->dict().lookupOrDefault
            (
                "implicitUSolve", implicitUCorrector_
            );
        momentumPredictor_ =
            solnControlPtr_->dict().lookupOrDefault
            (
                "momentumPredictor", false
            );
    }
    else if (correctorName == "PISOCorrector")
    {
        corr_ = corrector;
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        nonOrthCorr_ = corrector;
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::assembleUEqnLHS()
{
    volVectorField& U = U_();
    surfaceScalarField& phi = phi_();
    volScalarField& rho = rho_();
    rho = thermoPtr_->rho();

    fv::options& fvOptions = this->fvOptions();

    adjustSplitBCs(phi_(), U_(), U);

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

    UEqnLHS_.ref() += fvm::ddt(rho, U);

    UEqnLHS_->relax();

    fvOptions.constrain(UEqnLHS_.ref());
}


void Foam::fv::multiphaseSegregatedFlowSolver::assemblepEqn()
{
    volVectorField& U = U_();
    surfaceScalarField& phiv = phiv_();
    volScalarField& rho = rho_();
    volScalarField& p = *p_;
    const volScalarField& psi = thermoPtr_->psi();
    fv::options& fvOptions = this->fvOptions();

    if (frameAcceleration_.valid())
    {
        calculateghFields(buoyant_);
    }
    calculateFrameAccelerationContribution();

    if (!momentumPredictor_ && corr_ == 0 && nonOrthCorr_ == 0)
    {
        assembleUEqnLHS();
    }

    if (nonOrthCorr_ == 0)
    {
        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution
        p0_ = tmp<volScalarField>(new volScalarField(p));

        fvVectorMatrix& UEqnLHS(UEqnLHS_.ref());

        tmp<volScalarField> H1;
        if (solnControl().consistent())
        {
            tmp<volScalarField> A(UEqnLHS.A());
            // Stabilise the degenerate case where H1 ~= -A
            H1 = min(UEqnLHS.H1(), (1.0-0.01)*A());
            rAU_ = new volScalarField("(1|A(U))", 1/(A-H1()));
        }
        else
        {
            rAU_ = new volScalarField("(1|A(U))", 1/UEqnLHS.A());
        }

        adjustSplitBCs(phi_(), U_(), p);

        const volScalarField& rAU = rAU_();
        if (solnControl().consistent())
        {
            HbyA_ =
                new volVectorField
                (
                    constrainHbyA(rAU*(UEqnLHS.H() - H1*U), U, p)
                );
        }
        else
        {
            HbyA_ =
                new volVectorField
                (
                    constrainHbyA(rAU*UEqnLHS.H(), U, p)
                );
        }
        const volVectorField& HbyA = HbyA_();

        phivHbyA_ = new surfaceScalarField("phivHbyA", fvc::flux(HbyA));
        surfaceScalarField& phivHbyA(phivHbyA_.ref());
        rAUf_ = new surfaceScalarField(rAUfName_, fvc::interpolate(rAU));

        tmp<surfaceScalarField> phivg;
        if (buoyant_)
        {
            phivg =
                rAUf_()
               *(
                    this->faceBuoyancyForce(rAUfName_)
                  + thermoPtr_->surfaceTensionForce(U)
                )*mesh_.magSf();
        }
        else
        {
            phivg = rAUf_()*thermoPtr_->surfaceTensionForce(U)*mesh_.magSf();
        }

        phivHbyA += fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phiv, Uf_);
        fvOptions.makeRelative(phivHbyA);

        closedVolume_ = p.needReference();

        phivHbyA += phivg;

        // Update the fixedFluxPressure BCs to ensure flux consistency
        constrainPressure(p, U, phivHbyA, rAUf_(), fvOptions);

        pEqnComps_.setSize(thermoPtr_->thermos().size());

        forAll(thermoPtr_->thermos(), phasei)
        {
            const rhoThermo& thermoi =
                refCast<const rhoThermo>(thermoPtr_->thermos()[phasei]);
            tmp<volScalarField> trhoi = thermoi.rho();
            const volScalarField& rhoi = trhoi();

            if (!thermoi.isochoric())
            {
                pEqnComps_.set
                (
                    phasei,
                    new fvScalarMatrix
                    (
                        p, rho.dimensions()/dimTime*dimVolume
                    )
                );
                pEqnComps_[phasei] +=
                    fvc::ddt(rhoi)
                  + fvc::div(phiv, rhoi) - fvc::Sp(fvc::div(phiv), rhoi);
                // Compensate the div(meshPhi) in ddt(rhoi), since it cancels
                // in the two div terms above and the incompressible part of the
                // equation uses absolute fluxes
                if (mesh_.moving())
                {
                    pEqnComps_[phasei] -=
                        fvc::div
                        (
                            mesh_.phi(),
                            rhoi,
                            fvc::divSchemeName(phiv, rhoi)
                        );
                }

                if (!thermoi.incompressible())
                {
                    pEqnComps_[phasei] += thermoi.psi()*correction(fvm::ddt(p));
                }
            }
        }

        if (!thermoPtr_->incompressible())
        {
            compressibility_() = fvc::domainIntegrate(psi);
        }

        pDDtEqn_ =
            new fvScalarMatrix
            (
              - fvOptions(psi/rho, p)
              + fvc::div(phivHbyA)
            );
    }

    tmp<fv::laplacianScheme<scalar, scalar>> pLaplacianScheme
    (
        fv::laplacianScheme<scalar, scalar>::New
        (
            mesh_,
            p.db(),
            mesh_.schemes().laplacianScheme
            (
                "laplacian("+rAUfName_+","+p.name()+')'
            ),
            "grad("+p.name()+')'
        )
    );

    // Only apply the snGrad scheme's limiting to the combination
    // p-rho*gh rather than p directly.
    tmp<surfaceScalarField> limitedSnGradp;
    if (thermoPtr_->buoyant())
    {
        limitedSnGradp =
            rAUf_()
           *(
                pLaplacianScheme().snGradientScheme().snGrad(p-rho*gh_())
              + pLaplacianScheme().snGradientScheme().snGrad(rho_()*gh_())
            )*mesh_.magSf();
    }
    else
    {
        limitedSnGradp =
            rAUf_()
           *(pLaplacianScheme().snGradientScheme().snGrad(p))*mesh_.magSf();
    }

    // Solve pressure
    pEqnIncomp_ =
        new fvScalarMatrix
        (
            pDDtEqn_()
          - correction(pLaplacianScheme.ref().fvmLaplacian(rAUf_(), p))
          - fvc::div(limitedSnGradp())
        );

    // pEqnIncomp.flux() will have faces already masked by the GIB mask, so
    // we do this as well with limitedSnGradp to keep everything in balance
    fvc::applyFaceMaskTo(limitedSnGradp.ref());
    *pEqnIncomp_.ref().faceFluxCorrectionPtr() -= limitedSnGradp;

    if (finalIter_["nonOrthogonalCorrector:"+p_->name()])
    {
        rAUf_.clear();
    }

    pEqn_ = new fvScalarMatrix(pEqnIncomp_());

    forAll(thermoPtr_->thermos(), phasei)
    {
        if (pEqnComps_.set(phasei))
        {
            tmp<fvScalarMatrix> hmm
            (
                (
                    max(thermoPtr_->alphas()[phasei], scalar(0))
                   /thermoPtr_->thermos()[phasei].rho()
                )
               *pEqnComps_[phasei]
            );

            pEqn_.ref() += hmm;
        }
    }

    // Do not set reference in (transient), compressible case as there is a
    // dp/dt term
    if (p.needReference() && thermoPtr_->incompressible())
    {
        pEqn_->setReference
        (
            pressureControl_->refCell(),
            pressureControl_->refValue()
        );
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fv::multiphaseSegregatedFlowSolver::faceBuoyancyForce
(
    const word& pDiffName, bool includeSnGradP
) const
{
    // Overriding this function to cancel a term; this is not the full
    // buoyancy force but has a snGrad(rho*gh) removed. Will be refactored.

    tmp<volScalarField> tbRho = buoyantRho();
    const volScalarField& bRho = tbRho();

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
                "laplacian("+pDiffName+","+p().name()+')'
            ),
            "grad("+p().name()+')'
        )
    );

    tmp<surfaceScalarField> fbSource =
        - ghf_()*pLaplacianScheme().snGradientScheme().snGrad(bRho)
        + pLaplacianScheme().snGradientScheme().snGrad(bRho*gh_());

    // Note: if includeSnGradP==true, this leaves it unbalanced, thereby
    // including it
    if (includeSnGradP)
    {
        fbSource.ref() -=
            pLaplacianScheme().snGradientScheme().snGrad(p()-bRho*gh_())
          + pLaplacianScheme().snGradientScheme().snGrad(bRho*gh_());
    }
    correctInactiveGIBZoneFaces(fbSource.ref());

    return fbSource;
}


tmp<fvScalarMatrix>
Foam::fv::multiphaseSegregatedFlowSolver::assembleScalarMatrix
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
                finalIter_[solverObject::outerCorrectorName]
             && finalIter_["PISOCorrector"]
             && finalIter_["nonOrthogonalCorrector:"+p_->name()]
            );
        return pEqn_;
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }

}


tmp<fvVectorMatrix>
Foam::fv::multiphaseSegregatedFlowSolver::assembleVectorMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if
    (
        (momentumPredictor_ && fieldName == U_->name()+"Predictor")
     || fieldName == U_->name()
    )
    {
        if (fieldName == U_->name()+"Predictor")
        {
            assembleUEqnLHS();
        }

        tmp<surfaceScalarField> tsnGrad;
        if (buoyant_)
        {
            tsnGrad = -this->faceBuoyancyForce(rAUfName_, true);
        }
        else
        {
            tsnGrad = fvc::snGrad(*p_);
        }
        tsnGrad.ref() -= thermoPtr_->surfaceTensionForce(U_());

        if
        (
            fieldName == U_->name()+"Predictor"
         || (implicitUSolve_ && corr_ == 0)
         || (implicitUCorrector_ && corr_ != 0)
        )
        {
            return UEqnLHS_() + fvc::reconstruct(tsnGrad*mesh_.magSf());
        }
        else
        {
            tmp<fvVectorMatrix> UEqn
            (
                UEqnLHS_() + fvc::reconstruct(tsnGrad*mesh_.magSf())
            );
            U_() = UEqn().H()/UEqn().A();
            U_->correctBoundaryConditions();
            return tmp<fvVectorMatrix>();
        }
    }
    else
    {
        return tmp<fvVectorMatrix>();
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();
    if (solveName == p_->name())
    {
        if (verbose_)
        {
            Info<< "p min/max: "
                << min(*p_).value() << " "
                << max(*p_).value() << endl;
        }

        if (finalIter_["nonOrthogonalCorrector:"+p_->name()])
        {
            // Final non-orthogonal corrector

            volScalarField& p = *p_;
            volScalarField& rho = rho_();

            forAll(thermoPtr_->thermos(), phasei)
            {
                if (pEqnComps_.set(phasei))
                {
                    thermoPtr_->dgdts()[phasei] =
                        pos0(thermoPtr_->alphas()[phasei])
                       *(pEqnComps_[phasei] & p)
                       /thermoPtr_->thermos()[phasei].rho();
                }
                else
                {
                    thermoPtr_->dgdts()[phasei] *= 0;
                }
            }
            pEqnComps_.clear();

            // Calculate the conservative flux
            phiv_() = phivHbyA_ + pEqnIncomp_->flux();
            pEqnIncomp_.clear();
            if (!isStatic())
            {
                fvc::makeRelative(phiv_(), *U_);
            }

            pressureControl_->limit(p);

            if (closedVolume_)
            {
                if (thermoPtr_->incompressible())
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
                                thermoPtr_->rho() + thermoPtr_->psi()*(p-p0_())
                            )
                        )/compressibility_();
                }
            }

            if (!thermoPtr_->incompressible())
            {
                //TODO: make single phase thermo also work with dp instead of drho to be consistent with this
                thermoPtr_->correctRho(p - p0_());
            }

            rho = thermoPtr_->rho();
        }
        else
        {
            pEqnIncomp_.clear();
        }
    }
    else if
    (
        (momentumPredictor_ && solveName == U_->name()+"Predictor")
     || solveName == U_->name()
    )
    {
        // U corrector
        volVectorField& U = U_();

        fvOptions.correct(U_());

        if (verbose_)
        {
            Info<< "mag(U) min/max: "
                << min(mag(U_())).value() << " "
                << max(mag(U_())).value() << endl;
        }

        if (Uf_.valid() && solveName == U_->name())
        {
            Uf_() = fvc::interpolate(U);
            surfaceVectorField n(mesh_.Sf()/mesh_.magSf());
            Uf_() +=
                n
               *(
                    fvOptions.absolute(fvc::absolute(phiv_(), U))/mesh_.magSf()
                  - (n & Uf_())
                );
        }
    }
    else if (solveName == "correctPhi")
    {
        // Perform delayed initialisation (needed to be done after first mesh
        // move)
        if (phiNeedsInit_)
        {
            this->fvOptions().makeRelative(phiv_());
            tmp<surfaceScalarField> rhof(fvc::interpolate(rho_()));
            this->fvOptions().makeRelative(rhof, phi_());
            if (!isStatic())
            {
                fvc::makeRelative(phiv_(), U_());
                fvc::makeRelative(phi_(), rho_(), U_());
            }
            phiNeedsInit_ = false;
        }
        if (!isStatic() && !divU_.valid())
        {
            divU_ =
                fvc::div
                (
                    this->fvOptions().absolute(fvc::absolute(phiv_(), U_()))
                );
        }

        if (needsCorrectPhi_)
        {
            correctPhi(true);
            calculateghFields(buoyant_);
            calculateFrameAccelerationContribution();
            needsCorrectPhi_ = false;
        }
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::solveContinuity()
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


void Foam::fv::multiphaseSegregatedFlowSolver::printIncoContinuityErrors
(
    const word& regionName
)
{
    volScalarField contErr(fvc::div(phi_()));

    scalar sumLocalContErr = mesh_.time().deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    scalar globalContErr = mesh_.time().deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value();
    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::multiphaseSegregatedFlowSolver::printCompContinuityErrors
(
    const word& regionName
)
{
    const volScalarField& rho = rho_();

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr, globalContErr;

    sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rho - thermoPtr_->rho()))/totalMass
    ).value();

    globalContErr =
    (
        fvc::domainIntegrate(rho - thermoPtr_->rho())/totalMass
    ).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors (" << regionName << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::fv::multiphaseSegregatedFlowSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        if (!isStatic())
        {
            // Store divU from this mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phiv has the
            // same divergence after the mesh update
            divU_ =
                tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "divU",
                        fvc::div(fvc::absolute(phiv_(), U_()))
                    )
                );
        }

        if (finalIter)
        {
            if (!fv::localEulerDdt::enabled(mesh_))
            {
                Info<< endl;
                // Calculate Courant for next time step
                maxCo_ =
                    CourantNumber(mesh_, mesh_.time().deltaTValue(), geometricOneField(), phiv_());
            }
        }
    }
}


void Foam::fv::multiphaseSegregatedFlowSolver::correctPhi(bool defaultCorrectPhi)
{
    if
    (
        solnControlPtr_->dict().lookupOrDefault<Switch>
        (
            "correctPhi", defaultCorrectPhi
        )
    )
    {
        // Calculate absolute flux from the mapped surface velocity
        phiv_().forceAssign(mesh_.Sf() & fvc::interpolate(U_()));

        CorrectPhi
        (
            U_(),
            phiv_(),
            *p_,
            dimensionedScalar("rAUf", dimTime/rho_().dimensions(), 1),
            divU_(),
            *solnControlPtr_,
            // For correctPhi, always use a p reference (if boundary conditions
            // require it), but only set it on one processor
            getFirstCell(mesh_).first()
        );

        // Make the fluxes relative to the mesh-motion and MRF
        this->fvOptions().makeRelative(phiv_());
        fvc::makeRelative(phiv_(), U_());
    }
}


bool Foam::fv::multiphaseSegregatedFlowSolver::movePoints()
{
    if (mesh_.changing()) limitVolumes(this->dict(), mesh_);
    calculateghFields(buoyant_);
    needsCorrectPhi_ = true;
    calculateFrameAccelerationContribution();
    return true;
}


void Foam::fv::multiphaseSegregatedFlowSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (mesh_.changing()) limitVolumes(this->dict(), mesh_);
    needsCorrectPhi_ = true;
}


// ************************************************************************* //
