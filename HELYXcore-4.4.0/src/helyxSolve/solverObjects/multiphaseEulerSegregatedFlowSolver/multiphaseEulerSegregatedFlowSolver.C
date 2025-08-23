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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseEulerSegregatedFlowSolver.H"
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
#include "cfdTools/general/limitVolumes/limitVolumes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseEulerSegregatedFlowSolver, 0);
}
}

makeFvSolverOption(multiphaseEulerSegregatedFlowSolver);

const word Foam::fv::multiphaseEulerSegregatedFlowSolver::rAUfName_ = "rAUf";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseEulerSegregatedFlowSolver::multiphaseEulerSegregatedFlowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    flowSolver(name, obr, dict),
    thermoPtr_(nullptr),
    printContErr_(true),
    solnControlPtr_(nullptr),
    buoyant_(false)
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

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ = &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr));
    thermoPtr_->addPhasicVariable("U");

    buoyant_ = thermoPtr_->buoyant();
    if (thermoPtr_->distinctBuoyancy())
    {
        FatalErrorInFunction
            << "This solver does not support distinct buoyancy." << nl << endl;
    }

    if (buoyant_)
    {
        initializeGravityHref();
    }
    initializeFrameAcceleration(dict);
    calculateFrameAccelerationContribution();

    Info<< "Creating eulerianPhaseSystem\n" << endl;
    phaseSystem_.set
    (
        eulerianMultiphaseSystem::New(mesh_, *thermoPtr_).ptr()
    );

    // Create phi field for boundary conditions, postprocessing etc
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
            linearInterpolate(phaseSystem_().rho())*phaseSystem_().phiv()
        )
    );

    pressureControl_.set
    (
        new pressureControl
        (
            p,
            phaseSystem_().rho()(),
            solnControlPtr_->dict(),
            thermoPtr_->pRefValue(),
            p.needReference() && thermoPtr_->incompressible()
        )
    );

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

    Info<< "Creating averaged U field\n" << endl;
    U_.set(new volVectorField("U", phaseSystem_->U()));

    cumulativeContErr_ = 0;

    limitVolumes(this->dict(), mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseEulerSegregatedFlowSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    eulerianMultiphaseSystem::phaseModelList& phases = phaseSystem_().phases();

    // Solves
    solveNames = {p_->name()};

    // Dependencies
    optionalDependencies.insert(p_->name(), {"fvMesh"});

    DynamicList<word> UNames, heNames;
    forAll(phases, phasei)
    {
        eulerianPhaseModel& phase = phases[phasei];
        UNames.append(phase.U()().name());
    }

    forAll(phases, phasei)
    {
        solveNames.append(UNames[phasei]);
        requiredDependencies.insert(UNames[phasei], {p_->name()});
    }

    requiredDependencies.insert(p_->name(), {"alphas"});

    optionalDependencies.insert("turbulence", solveNames);

    // Correctors

    // Add all solves to outer corrector
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);

    // PISO (inner) correctors
    DynamicList<word> PISOMembers(UNames);
    PISOMembers.append(p_->name());
    correctorMembers.insert("PISOCorrector", PISOMembers);

    // Non-orthogonal correctors
    correctorMembers.insert
    (
        "nonOrthogonalCorrector:" + p_->name(), {p_->name()}
    );
}


bool Foam::fv::multiphaseEulerSegregatedFlowSolver::initialise()
{
    // Calculate initial Courant number
    maxCo_ =
        CourantNumber
        (
            mesh_,
            mesh_.time().deltaTValue(),
            phaseSystem_().rho()(),
            phaseSystem_().phiv()
        );

    return true;
}


scalar Foam::fv::multiphaseEulerSegregatedFlowSolver::getMaxTimeStep()
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


bool Foam::fv::multiphaseEulerSegregatedFlowSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        return false;
    }
    else if (correctorName == solverObject::outerCorrectorName)
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


void Foam::fv::multiphaseEulerSegregatedFlowSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == "PISOCorrector")
    {
        corr_ = corrector;
    }
    else if (correctorName == "nonOrthogonalCorrector:"+p_->name())
    {
        nonOrthCorr_ = corrector;
    }
}


void Foam::fv::multiphaseEulerSegregatedFlowSolver::assembleUEqnsLHSs()
{
    autoPtr<eulerianPhaseSystem::momentumTransferTable>
        momentumTransferPtr(phaseSystem_().momentumTransfer());

    eulerianPhaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    UEqns_.clear();
    UEqns_.setSize(phaseSystem_().phases().size());

    forAll(phaseSystem_().phases(), phasei)
    {
        eulerianPhaseModel& phase = phaseSystem_().phases()[phasei];
        phase.alphaPhi() = fvc::interpolate(phase.rho())*phase.alphaPhiv();
    }
    phaseSystem_().correct();

    forAll(phaseSystem_().phases(), phasei)
    {
        eulerianPhaseModel& phase = phaseSystem_().phases()[phasei];

        const volScalarField& alpha = phase;
        tmp<volScalarField> trho(phase.rho());
        const volScalarField& rho = trho();
        volVectorField& U = phase.U();

        UEqns_.set
        (
            phasei,
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
                *momentumTransfer[phase.name()]
              + fvOptions()(alpha, rho, U)
            )
        );

        UEqns_[phasei].relax();
        fvOptions().constrain(UEqns_[phasei]);
        fvOptions().correct(U);
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::fv::multiphaseEulerSegregatedFlowSolver::assemblepEqn()
{
    eulerianPhaseSystem::phaseModelList& phases = phaseSystem_().phases();
    volScalarField& p(*p_);
    surfaceScalarField& phiv = phaseSystem_().phiv();

    if (frameAcceleration_.valid())
    {
        calculateghFields(buoyant_);
    }
    calculateFrameAccelerationContribution();

    if (corr_ == 0 && nonOrthCorr_ == 0)
    {
        assembleUEqnsLHSs();

        alphafs_.clear();
        rAUs_.clear();
        alpharAUfs_.clear();
        alphafs_.setSize(phases.size());
        rAUs_.setSize(phases.size());
        alpharAUfs_.setSize(phases.size());

        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];
            const volScalarField& alpha = phase;

            alphafs_.set(phasei, fvc::interpolate(alpha).ptr());
            alphafs_[phasei].rename("pEqn" + alphafs_[phasei].name());

            rAUs_.set
            (
                phasei,
                new volScalarField
                (
                    IOobject::groupName("rAU", phase.name()),
                    1.0
                   /(
                        UEqns_[phasei].A()
                      + max(phase.residualAlpha() - alpha, scalar(0))
                       *phase.rho()/mesh_.time().deltaT()
                    )
                )
            );

            alpharAUfs_.set
            (
                phasei,
                (
                    fvc::interpolate
                    (
                        max(alpha, phase.residualAlpha())*rAUs_[phasei]
                    )
                ).ptr()
            );
        }

        // Lift, wall-lubrication and turbulent diffusion fluxes
        phiFs_.clear();
        phiFs_.setSize(phases.size());
        {
            autoPtr<PtrList<volVectorField>> Fs = phaseSystem_().Fs();

            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];

                if (Fs().set(phasei))
                {
                    phiFs_.set
                    (
                        phasei,
                        new surfaceScalarField
                        (
                            IOobject::groupName("phiF", phase.name()),
                            fvc::flux(rAUs_[phasei]*Fs()[phasei])
                        )
                    );
                }
            }
        }
        {
            autoPtr<PtrList<surfaceScalarField>> phiDs =
                phaseSystem_().phiDs(rAUs_);

            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];

                if (phiDs().set(phasei))
                {
                    if (phiFs_.set(phasei))
                    {
                        phiFs_[phasei] += phiDs()[phasei];
                    }
                    else
                    {
                        phiFs_.set
                        (
                            phasei,
                            new surfaceScalarField
                            (
                                IOobject::groupName("phiF", phase.name()),
                                phiDs()[phasei]
                            )
                        );
                    }
                }
            }
        }
    }

    if (nonOrthCorr_ == 0)
    {
        volScalarField rho("rho", phaseSystem_().rho());

        HbyAs_.clear();
        HbyAs_.setSize(phases.size());

        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];
            const volScalarField& alpha = phase;

            // Correct fixed-flux BCs to be consistent with the velocity BCs
            fvOptions().correctBoundaryFlux(phase.U(), phase.phiv());

            HbyAs_.set
            (
                phasei,
                new volVectorField
                (
                    IOobject::groupName("HbyA", phase.name()),
                    phase.U()
                )
            );

            HbyAs_[phasei] =
                rAUs_[phasei]
               *(
                    UEqns_[phasei].H()
                  + max(phase.residualAlpha() - alpha, scalar(0))
                   *phase.rho()*phase.U().oldTime()/mesh_.time().deltaT()
                );
        }

        phigFs_.clear();
        phigFs_.setSize(phases.size());
        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];

            if (buoyant_)
            {
                phigFs_.set
                (
                    phasei,
                    (
                        -alpharAUfs_[phasei]
                        *(
                            (
                                this->faceBuoyancyForce(rAUfName_)
                              + thermoPtr_->surfaceTensionForce
                                (
                                    phase.U(), phase.volFrac()
                                )
                            )*mesh_.magSf()
                          + (fvc::interpolate(phase.rho() - rho))
                           *(g_() & mesh_.Sf())
                        )
                    ).ptr()
                );
            }
            else
            {
                phigFs_.set
                (
                    phasei,
                    (
                        -alpharAUfs_[phasei]
                        *thermoPtr_->surfaceTensionForce
                         (
                            phase.U(), phase.volFrac()
                         )*mesh_.magSf()
                    ).ptr()
                );
            }

            if (phiFs_.set(phasei))
            {
                phigFs_[phasei] += phiFs_[phasei];
            }
        }

        phivHbyAs_.clear();
        phivHbyAs_.setSize(phases.size());

        phivHbyA_ =
            new surfaceScalarField
            (
                IOobject
                (
                    "phiHbyA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimArea*dimVelocity, 0)
            );

        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];
            const volScalarField& alpha = phase;

            // ddtPhiCorr filter -- only apply in pure(ish) phases
            surfaceScalarField alphafBar
            (
                fvc::interpolate(fvc::average(alphafs_[phasei]))
            );
            surfaceScalarField phiCorrCoeff(pos0(alphafBar - 0.99));

            surfaceScalarField::Boundary& phiCorrCoeffBf =
                phiCorrCoeff.boundaryFieldRef();

            forAll(mesh_.boundary(), patchi)
            {
                // Set ddtPhiCorr to 0 on non-coupled boundaries
                if
                (
                    !mesh_.boundary()[patchi].coupled()
                 || isA<cyclicAMIFvPatch>(mesh_.boundary()[patchi])
                )
                {
                    phiCorrCoeffBf[patchi] = 0;
                }
            }

            phivHbyAs_.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject::groupName("phiHbyA", phase.name()),
                    fvc::flux(HbyAs_[phasei])
                  + phiCorrCoeff
                   *fvc::interpolate
                    (
                        alpha.oldTime()*phase.rho()().oldTime()*rAUs_[phasei]
                    )
                   *(
                        fvOptions().absolute(phase.phiv().oldTime())
                      - fvc::flux(phase.U().oldTime())
                    )/mesh_.time().deltaT()
                  - phigFs_[phasei]
                )
            );

            forAllConstIter(eulerianPhaseSystem::KdTable, phaseSystem_().Kds(), KdIter)
            {
                const volScalarField& K(*KdIter());

                const eulerianPhasePair& pair
                (
                    phaseSystem_().phasePairs()[KdIter.key()]
                );

                const eulerianPhaseModel* phase1 = &pair.phase1();
                const eulerianPhaseModel* phase2 = &pair.phase2();

                forAllConstIter(eulerianPhasePair, pair, iter)
                {
                    if (phase1 == &phase)
                    {
                        phivHbyAs_[phasei] +=
                            fvc::interpolate(rAUs_[phasei]*K)
                           *fvOptions().absolute(phase2->phiv());

                        HbyAs_[phasei] += rAUs_[phasei]*K*phase2->U();
                    }

                    Swap(phase1, phase2);
                }
            }

            phivHbyA_.ref() += alphafs_[phasei]*phivHbyAs_[phasei];
        }

        fvOptions().makeRelative(phivHbyA_.ref());

        // Construct pressure "diffusivity"
        rAUf_ =
            new surfaceScalarField
            (
                IOobject
                (
                    rAUfName_,
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(rAUfName_, dimensionSet(-1, 3, 1, 0, 0), 0)
            );
        surfaceScalarField& rAUf = rAUf_.ref();

        forAll(phases, phasei)
        {
            rAUf += alphafs_[phasei]*alpharAUfs_[phasei];
        }
        rAUf = mag(rAUf);

        // Update the fixedFluxPressure BCs to ensure flux consistency
        {
            surfaceScalarField::Boundary phivb(phiv.boundaryField());
            phivb = 0;
            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];
                phivb +=
                    alphafs_[phasei].boundaryField()*phase.phiv().boundaryField();
            }

            setSnGrad<fixedFluxPressureFvPatchScalarField>
            (
                p.boundaryFieldRef(),
                (
                    phivHbyA_().boundaryField() - phivb
                )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
            );
        }

        pEqnComps_.clear();
        pEqnComps_.setSize(phases.size());
        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];
            const volScalarField& alpha = phase;
            volScalarField& rho = phase.thermo().rho();

            if (phase.compressible())
            {
                if (solnControlPtr_->transonic())
                {
                    surfaceScalarField phid
                    (
                        IOobject::groupName("phid", phase.name()),
                        fvc::interpolate(phase.thermo().psi())*phase.phiv()
                    );

                    pEqnComps_.set
                    (
                        phasei,
                        (
                            (
                                phase.continuityError()
                              - fvc::Sp
                                (
                                    fvc::ddt(alpha)
                                  + fvc::div(phase.alphaPhiv()),
                                    rho
                                )
                            )/rho
                          + correction
                            (
                                (alpha/rho)
                               *(
                                    phase.thermo().psi()*fvm::ddt(p)
                                  + fvm::div(phid, p)
                                  - fvm::Sp(fvc::div(phid), p)
                                )
                            )
                        ).ptr()
                    );

                    deleteDemandDrivenData
                    (
                        pEqnComps_[phasei].faceFluxCorrectionPtr()
                    );
                    pEqnComps_[phasei].relax();
                }
                else
                {
                    pEqnComps_.set
                    (
                        phasei,
                        (
                            (
                                phase.continuityError()
                              - fvc::Sp
                                (
                                    (
                                        fvc::ddt(alpha)
                                      + fvc::div(phase.alphaPhiv())
                                    ),
                                    rho
                                )
                            )/rho
                          + (alpha*phase.thermo().psi()/rho)
                           *correction(fvm::ddt(p))
                        ).ptr()
                    );
                }
            }
            else
            {
                pEqnComps_.set
                (
                    phasei,
                    fvm::Su(-(fvOptions()(alpha, rho)&rho)/rho, p).ptr()
                );
            }

            if (phaseSystem_().transfersMass(phase))
            {
                if (pEqnComps_.set(phasei))
                {
                    pEqnComps_[phasei] -= phaseSystem_().dmdt(phase)/rho;
                }
                else
                {
                    pEqnComps_.set
                    (
                        phasei,
                        fvm::Su(-phaseSystem_().dmdt(phase)/rho, p)
                    );
                }
            }
        }

        // Cache p prior to solve for density update
        p0_ = new volScalarField(p);
    }

    // Construct the transport part of the pressure equation
    pEqnIncomp_ =
        fvc::div(phivHbyA_())
      - fvm::laplacian(rAUf_(), p);

    tmp<fvScalarMatrix> pEqn(new fvScalarMatrix(pEqnIncomp_()));

    forAll(phases, phasei)
    {
        if (pEqnComps_.set(phasei))
        {
            pEqn.ref() += pEqnComps_[phasei];
        }
    }

    return pEqn;
}


tmp<fvScalarMatrix>
Foam::fv::multiphaseEulerSegregatedFlowSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if (fieldName == p_->name())
    {
        finalSolve =
            (
                finalIter_[solverObject::outerCorrectorName]
             && finalIter_["PISOCorrector"]
             && finalIter_["nonOrthogonalCorrector:"+p_->name()]
            );
        return assemblepEqn();
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }
}


void Foam::fv::multiphaseEulerSegregatedFlowSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();

    if (solveName == p_->name())
    {
        volScalarField& p(*p_);
        if (verbose_)
        {
            Info<< "p min/max: "
                << min(p).value() << " "
                << max(p).value() << endl;
        }

        if (finalIter_["nonOrthogonalCorrector:"+p_->name()])
        {
            surfaceScalarField& phiv = phaseSystem_().phiv();
            eulerianPhaseSystem::phaseModelList& phases = phaseSystem_().phases();

            // Correct fluxes and velocities on last non-orthogonal iteration
            phiv = phivHbyA_() + pEqnIncomp_().flux();

            mSfGradp_ =
                new surfaceScalarField("mSfGradp", pEqnIncomp_().flux()/rAUf_());

            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];

                phase.phiv() =
                    phivHbyAs_[phasei] + alpharAUfs_[phasei]*mSfGradp_();

                // Set the phase dilatation rates
                if (pEqnComps_.set(phasei))
                {
                    phase.divU(-pEqnComps_[phasei] & p);
                    volScalarField& dgdt = thermoPtr_->dgdts()[phasei];
                    volScalarField& alpha = phases[phasei].volFrac();
                    dgdt = -phase.divU()()/max(alpha, scalar(1e-4));
                }
            }
            phivHbyAs_.clear();

            // Optionally relax pressure for velocity correction
            p.relax();

            mSfGradp_ = pEqnIncomp_().flux()/rAUf_;
            pEqnIncomp_.clear();

            // Update and limit the static pressure
            pressureControl_().limit(p);

            // Update densities from change in p
            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];
                phase.thermo().rho() += phase.thermo().psi()*(p - p0_());
            }
            p0_.clear();

        }
    }
    else
    {
        const auto i = solveName.rfind('.');
        word solveMember, solveGroup;
        if (i == std::string::npos || i == 0)
        {
            solveMember = solveName;
            solveGroup = word::null;
        }
        else
        {
            solveMember = solveName.substr(0, i);
            solveGroup = solveName.substr(i+1);
        }

        if (solveMember == "U")
        {
            // U corrector
            eulerianPhaseSystem::phaseModelList& phases = phaseSystem_().phases();
            eulerianPhaseModel& phase = phases[solveGroup];
            const label phasei = phase.index();

            phase.U() =
                HbyAs_[phasei]
              + fvc::reconstruct
                (
                    alpharAUfs_[phasei]*mSfGradp_()
                  - phigFs_[phasei]
                );
            phase.U().correctBoundaryConditions();
            fvOptions.correct(phase.U());
        }
    }
}


void Foam::fv::multiphaseEulerSegregatedFlowSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == "PISOCorrector")
    {
        HbyAs_.clear();
        phigFs_.clear();
        mSfGradp_.clear();
        if (finalIter)
        {
            UEqns_.clear();
            alpharAUfs_.clear();
        }

        // Update averaged U for post-processing etc.
        U_() = phaseSystem_().U();
    }
    else if (correctorName == solverObject::outerCorrectorName)
    {
        if (finalIter)
        {
            Info<< endl;
            // Calculate Courant for next time step
            maxCo_ =
                CourantNumber
                (
                    mesh_,
                    mesh_.time().deltaTValue(),
                    geometricOneField(),
                    phaseSystem_().phiv()
                );
        }
    }
}


bool Foam::fv::multiphaseEulerSegregatedFlowSolver::movePoints()
{
    FatalErrorInFunction
        << "This solver does not support dynamic meshes."
        << exit(FatalError);
    calculateghFields(buoyant_);
    calculateFrameAccelerationContribution();
    return true;
}


void Foam::fv::multiphaseEulerSegregatedFlowSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    FatalErrorInFunction
        << "This solver does not support dynamic meshes."
        << exit(FatalError);
}


// ************************************************************************* //
