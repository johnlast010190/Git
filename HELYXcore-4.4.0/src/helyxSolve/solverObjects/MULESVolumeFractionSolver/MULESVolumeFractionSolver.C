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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "MULESVolumeFractionSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "fvMatrices/solvers/MULES/CMULES.H"
#include "algorithms/subCycle/subCycle.H"
#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"
#include "eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "eulerianMultiphaseSystem/eulerianMultiphaseSystem.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(MULESVolumeFractionSolver, 0);
}
}

makeFvSolverOption(MULESVolumeFractionSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MULESVolumeFractionSolver::MULESVolumeFractionSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    thermoPtr_(nullptr),
    phaseSystemPtr_(nullptr),
    passiveIndex_(-1),
    phivPtr_(nullptr),
    phiPtr_(nullptr),
    maxAlphaCo_(0),
    solnControlPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::MULESVolumeFractionSolver::read(const dictionary& dict)
{
    if (dict.found("interfaceCompression"))
    {
        cAlphas_.reset(new HashTable<scalar>());
        const dictionary& tp(dict.subDict("interfaceCompression"));
        for (const word& k : tp.toc())
        {
            cAlphas_->insert(k, readScalar(tp.lookup(k)));
        }
    }
}


bool Foam::fv::MULESVolumeFractionSolver::initialise()
{
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    // Check for eulerianPhaseSystem (present in Euler-Euler case)
    phaseSystemPtr_ =
        obr_.lookupObjectRefPtr<eulerianMultiphaseSystem>
        (
            eulerianPhaseSystem::typeName
        );

    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr_));

    passiveIndex_ = thermoPtr_->fractions().passiveIndex();

    read(dictionary(dict()));

    phivPtr_ = &mesh_.lookupObject<surfaceScalarField>("phiv");

    // U is not used in Euler-Euler case; phi is only created for
    // postprocessing, so we treat it as optional from the point of view of
    // this solver object
    if (phaseSystemPtr_)
    {
        phiPtr_ = mesh_.lookupObjectRefPtr<surfaceScalarField>("phi");
        UPtr_ = nullptr;
    }
    else
    {
        phiPtr_ = &mesh_.lookupObjectRef<surfaceScalarField>("phi");
        UPtr_ = &mesh_.lookupObject<volVectorField>("U");
    }

    updateTimeSchemeInfo();

    return true;
}


void Foam::fv::MULESVolumeFractionSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    const PtrList<volScalarField>& alphas = thermoPtr_->alphas();
    DynamicList<word> alphaNames;
    forAll(alphas, i)
    {
        alphaNames.append(alphas[i].name());
    }
    solveNames.append("alphas");
    derivedFields.insert("alphas", alphaNames);
    optionalDependencies.insert("alphas", {"fvMesh"});
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


scalar Foam::fv::MULESVolumeFractionSolver::getMaxTimeStep()
{
    if (!LTS_)
    {
        // Try region-specific maxCo, otherwise revert to global
        const dictionary& maxCoDict
        (
            solnControlPtr_->dict().found("maxAlphaCo")
          ? solnControlPtr_->dict()
          : mesh_.time().controlDict()
        );
        if (maxCoDict.found("maxAlphaCo"))
        {
            maxAlphaCoDataPtr_ = Function1<scalar>::New("maxAlphaCo", maxCoDict);

            scalar maxMaxAlphaCo =
                maxAlphaCoDataPtr_->value(mesh_.time().timeIndex());
            return
                maxMaxAlphaCo/stabilise(maxAlphaCo_, SMALL)
               *mesh_.time().deltaTValue();
        }
        else
        {
            return GREAT;
        }
    }
    else
    {
        return GREAT;
    }
}


bool Foam::fv::MULESVolumeFractionSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    return true;
}


void Foam::fv::MULESVolumeFractionSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        const dictionary* alphaControls = nullptr;
        if (mesh().solution().found("alpha"))
        {
            alphaControls = &mesh().solution().solverDict("alpha");
        }
        // This can be moved to the read() function once deprecation is done
        if
        (
            alphaControls
         && isDeprecatedAlphaControl(*alphaControls, "nAlphaSubCycles")
        )
        {
            nAlphaSubCycles_ = alphaControls->lookup<label>("nAlphaSubCycles");
        }
        else
        {
            nAlphaSubCycles_ = dict().lookup<label>("nAlphaSubCycles");
        }
        if
        (
            alphaControls
         && isDeprecatedAlphaControl(*alphaControls, "nAlphaCorr")
        )
        {
            nAlphaCorr_ =
                alphaControls->lookupOrDefault<label>("nAlphaCorr", 1);
        }
        else
        {
            nAlphaCorr_ = dict().lookupOrDefault<label>("nAlphaCorr", 1);
        }
        if
        (
            alphaControls
         && isDeprecatedAlphaControl(*alphaControls, "MULESCorr")
        )
        {
            MULESCorr_ =
                alphaControls->lookupOrDefault<Switch>("MULESCorr", false);
        }
        else
        {
            MULESCorr_ = dict().lookupOrDefault<Switch>("MULESCorr", false);
        }
        if
        (
            alphaControls
         && isDeprecatedAlphaControl(*alphaControls, "MULESCorr")
        )
        {
            alphaApplyPrevCorr_ =
                alphaControls->lookupOrDefault<Switch>
                (
                    "alphaApplyPrevCorr", false
                );
        }
        else
        {
            alphaApplyPrevCorr_ =
                dict().lookupOrDefault<Switch>("alphaApplyPrevCorr", false);
        }
        if (alphaControls && isDeprecatedAlphaControl(*alphaControls, "cAlpha"))
        {
            cAlphaGlobal_ =
                alphaControls->lookupOrDefault<scalar>("cAlpha", 0.0);
        }
        else
        {
            cAlphaGlobal_ = dict().lookupOrDefault<scalar>("cAlpha", 0.0);
        }
        correctInflowOutflow_ =
            dict().lookupOrDefault<Switch>("correctInflowOutflow", true);

        if (!alphaApplyPrevCorr_)
        {
            alphaPhivCorr0_.clear();
        }

        updateTimeSchemeInfo();

        if (MULESCorr_)
        {
            forAll(thermoPtr_->alphas(), phasei)
            {
                if (phasei != passiveIndex_ || thermoPtr_->alphas().size() > 2)
                {
                    mesh().schemes().setFluxRequired
                    (
                        thermoPtr_->alphas()[phasei].name()
                    );
                }
            }
        }

    }
}


void Foam::fv::MULESVolumeFractionSolver::updateTimeSchemeInfo()
{
    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh_,
                mesh_.schemes().ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            crankNicolson_ = false;
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            crankNicolson_ = true;
            if (mesh_.time().subCycling())
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            ocCoeff =
                refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                .ocCoeff();
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Keep the alphaPhivs around
    if (!alphaPhivs_.size())
    {
        alphaPhivs_.setSize(thermoPtr_->phases().size());
        forAll(alphaPhivs_, phasei)
        {
            alphaPhivs_.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaPhiv",
                            thermoPtr_->phases()[phasei]
                        ),
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    (*phivPtr_)
                    *fvc::interpolate(thermoPtr_->alphas()[phasei])
                )
            );
        }
    }

    // Set the time blending factor, 1 for Euler
    cnCoeff_ = 1.0/(1.0 + ocCoeff);

    LTS_ = fv::localEulerDdt::enabled(mesh_);
}

void Foam::fv::MULESVolumeFractionSolver::calcSuSp
(
    const label phasei,
    autoPtr<volScalarField::Internal>& pSu,
    autoPtr<volScalarField::Internal>& pSp
)
{
    volScalarField& alpha = thermoPtr_->alphas()[phasei];
    const volScalarField& dgdt = thermoPtr_->dgdts()[phasei];

    pSu.set
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimVelocity/dimLength
        )
    );

    pSp.set
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dgdt.dimensions()
        )
    );

    volScalarField::Internal& Su(pSu());
    volScalarField::Internal& Sp(pSp());

    forAll(dgdt, celli)
    {
        if (dgdt[celli] < 0.0 && alpha[celli] > 0.0)
        {
            Sp[celli] =  dgdt[celli]*alpha[celli];
            Su[celli] = -dgdt[celli]*alpha[celli];
        }
        else if (dgdt[celli] > 0.0 && alpha[celli] < 1.0)
        {
            Sp[celli] = -dgdt[celli]*(1.0 - alpha[celli]);
            Su[celli] = 0.0;
        }
        else
        {
            Sp[celli] = 0.0;
            Su[celli] = 0.0;
        }
    }

    forAll(thermoPtr_->phases(), phase2i)
    {
        const volScalarField& alpha2 = thermoPtr_->alphas()[phase2i];

        if (&alpha2 != &alpha)
        {
            const scalarField& dgdt2 = thermoPtr_->dgdts()[phase2i];

            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*alpha2[celli];
                }
            }
        }
    }

    tmp<fvScalarMatrix> fvOpt = fvOptions()(alpha);
    if (fvOpt().hasLower() || fvOpt().hasUpper())
    {
        FatalErrorInFunction
            << "Only diagonal source terms are accepted in the fvOption "
            << "source for the alpha equations." << exit(FatalError);
    }
    Su.field() -= fvOpt().source()/mesh_.V();
    if (fvOpt().hasDiag())
    {
        Sp.field() += fvOpt().diag()/mesh_.V();
    }
}


void Foam::fv::MULESVolumeFractionSolver::solveAlphas
(
    const volScalarField::Internal& divU
)
{
    const surfaceScalarField& phiv(*phivPtr_);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    tmp<surfaceScalarField> phivCN(phiv);
    if (crankNicolson_)
    {
        phivCN = cnCoeff_*phiv + (1.0 - cnCoeff_)*phiv.oldTime();
        phivCN->rename(phiv.name());
    }

    const wordList& phases = thermoPtr_->phases();

    if (MULESCorr_)
    {
        volScalarField sumAlpha
        (
            "sumAlpha",
            phiv.db(),
            mesh_,
            dimless,
            0
        );
        surfaceScalarField sumAlphaPhiv
        (
            "sumAlphaPhiv",
            phiv.db(),
            mesh_,
            phiv.dimensions(),
            0
        );

        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_)
            {
                continue;
            }

            volScalarField& alpha = thermoPtr_->alphas()[phasei];

            const dictionary* MULEScontrols = &dict();
            if (mesh().solution().found(alpha.name()))
            {
                const dictionary& solnDict =
                    mesh().solution().solverDict(alpha.name());
                if
                (
                    isDeprecatedAlphaControl(solnDict, "nLimiterIter")
                 || isDeprecatedAlphaControl(solnDict, "smoothLimiter")
                 || isDeprecatedAlphaControl(solnDict, "extremaCoeff")
                )
                {
                    MULEScontrols = &solnDict;
                }
            }

            autoPtr<volScalarField::Internal> Su, Sp;
            calcSuSp(phasei, Su, Sp);

            fvScalarMatrix alphaEqn
            (
                (
                    LTS_
                  ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                  : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha)
                )
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh_,
                    phivCN,
                    upwind<scalar>(mesh_, phivCN)
                ).fvmDiv(phivCN, alpha)
             ==
                Su() + fvm::Sp(Sp() + divU, alpha)
            );

            alphaEqn.solve();

            tmp<surfaceScalarField> talphaPhivUD(alphaEqn.flux());
            alphaPhivs_[phasei] = talphaPhivUD();

            if
            (
                alphaApplyPrevCorr_
             && alphaPhivCorr0_.size()
             && alphaPhivCorr0_.set(phasei)
            )
            {
                MULES::correct
                (
                    alpha,
                    alphaPhivs_[phasei],
                    alphaPhivCorr0_[phasei],
                    1,
                    0,
                    MULEScontrols
                );

                alphaPhivs_[phasei] += alphaPhivCorr0_[phasei];
            }

            if (alphaApplyPrevCorr_)
            {
                // Cache the upwind-flux
                if (!alphaPhivCorr0_.size())
                {
                    alphaPhivCorr0_.setSize(thermoPtr_->phases().size());
                }
                // Note that this is updated after the correct step
                alphaPhivCorr0_.set(phasei, talphaPhivUD.ptr());
            }

            if (passiveIndex_ != -1)
            {
                sumAlpha += alpha;
                sumAlphaPhiv += alphaPhivs_[phasei];
            }
        }

        if (passiveIndex_ != -1)
        {
            volScalarField& alpha = thermoPtr_->alphas()[passiveIndex_];
            alpha.forceAssign(1-sumAlpha);
            alphaPhivs_[passiveIndex_] = phiv-sumAlphaPhiv;
        }
    }

    PtrList<surfaceScalarField> alphaPhivCorrs(phases.size());
    PtrList<surfaceScalarField> alphaPhivUns(phases.size());

    PtrList<volScalarField> alpha0(phases.size());

    for (label aCorr = 0; aCorr < nAlphaCorr_; aCorr++)
    {
        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_ && phases.size() <= 2)
            {
                continue;
            }

            const volScalarField& alpha = thermoPtr_->alphas()[phasei];

            // Use a phase-specific scheme only if supplied
            word scheme1("div(phiv," + alpha.name() + ')');
            word scheme2("div(phiv," + alpha.member() + ')');
            const word& alphaScheme
            (
               !mesh().schemes().divSchemes().found(scheme1)
             && mesh().schemes().divSchemes().found(scheme2)
             ? scheme2
             : scheme1
            );

            const dictionary* MULEScontrols = &dict();
            if (mesh().solution().found(alpha.name()))
            {
                const dictionary& solnDict =
                    mesh().solution().solverDict(alpha.name());
                if
                (
                    isDeprecatedAlphaControl(solnDict, "nLimiterIter")
                 || isDeprecatedAlphaControl(solnDict, "smoothLimiter")
                 || isDeprecatedAlphaControl(solnDict, "extremaCoeff")
                )
                {
                    MULEScontrols = &solnDict;
                }
            }

            tmp<surfaceScalarField> alphaPhivUn
            (
                new surfaceScalarField
                (
                    IOobject::groupName(phiv.name(), alpha.group()),
                    fvc::flux
                    (
                        phivCN(),
                        crankNicolson_
                      ? cnCoeff_*alpha + (1.0 - cnCoeff_)*alpha.oldTime()
                      : tmp<volScalarField>(alpha),
                        alphaScheme
                    )
                )
            );

            forAll(phases, phase2i)
            {
                const volScalarField& alpha2 = thermoPtr_->alphas()[phase2i];

                if (&alpha2 != &alpha)
                {
                    tmp<surfaceScalarField> phivr;
                    if (phaseSystemPtr_)
                    {
                        phivr =
                            phaseSystemPtr_->phases()[phasei].phiv()
                          - phaseSystemPtr_->phases()[phase2i].phiv();
                    }

                    bool cAlphaFound = false;
                    scalar cAlpha(0);
                    if (cAlphas_.valid())
                    {
                        auto cAlphaIter1
                        (
                            cAlphas_().find
                            (
                                phasePairKey
                                (
                                    phases[phasei],
                                    phases[phase2i],
                                    false
                                ).name()
                            )
                        );
                        auto cAlphaIter2
                        (
                            cAlphas_().find
                            (
                                phasePairKey
                                (
                                    phases[phase2i],
                                    phases[phasei],
                                    false
                                ).name()
                            )
                        );
                        cAlphaFound =
                            (
                                cAlphaIter1 != cAlphas_().end()
                             || cAlphaIter2 != cAlphas_().end()
                            );
                        if (cAlphaFound)
                        {
                            cAlpha =
                                (cAlphaIter1 != cAlphas_().end())
                              ? cAlphaIter1() : cAlphaIter2();
                        }
                    }
                    if (!cAlphaFound)
                    {
                        // Revert to global one if no pair-specific one found
                        cAlphaFound = (cAlphaGlobal_ > 0);
                        cAlpha = cAlphaGlobal_;
                    }

                    if (cAlphaFound)
                    {
                        tmp<surfaceScalarField> phivc;
                        if (phivr.valid())
                        {
                            phivc = (mag(phiv) + mag(phivr()))/mesh_.magSf();
                        }
                        else
                        {
                            phivc = mag(phiv)/mesh_.magSf();
                        }

                        tmp<surfaceScalarField> phivr2 =
                            min(cAlpha*phivc(), max(phivc()))
                           *thermoPtr_->nHatf(alpha, alpha2);
                        if (phivr.valid())
                        {
                            phivr.ref() += phivr2;
                        }
                        else
                        {
                            phivr = phivr2;
                        }
                    }

                    if (!phivr.valid())
                    {
                        // No compression terms from any source for this pair
                        continue;
                    }

                    // Use a phase-specific scheme only if supplied
                    word scheme1
                    (
                        "div(phivr,"+alpha2.name()+',' +alpha.name()+')'
                    );
                    word scheme2
                    (
                        "div(phivr,"+alpha2.member()+','+alpha.member()+')'
                    );
                    word oldAlpharName("div(phivrb,"+alpha.member()+')');
                    if (mesh().schemes().divSchemes().found(oldAlpharName))
                    {
                        static bool warned = false;
                        if (!warned)
                        {
                            DeprecationIOWarningInFunction
                            (
                                mesh().schemes().divSchemes(),
                                oldAlpharName,
                                "scheme name",
                                40300,
                                "Please use " + scheme2 + " instead."
                            );
                            warned = true;
                        }
                        if (!mesh().schemes().divSchemes().found(scheme2))
                        {
                            scheme2 = oldAlpharName;
                        }
                    }
                    const word& alpharScheme
                    (
                       !mesh().schemes().divSchemes().found(scheme1)
                     && mesh().schemes().divSchemes().found(scheme2)
                      ? scheme2
                      : scheme1
                    );
                    alphaPhivUn.ref() +=
                        fvc::flux
                        (
                            -fvc::flux(-phivr, alpha2, alpharScheme),
                            alpha,
                            alpharScheme
                        );
                }
            }

            // Remove interface compression at non-coupled boundaries
            if (phaseSystemPtr_)
            {
                correctInflowOutflow
                (
                    alphaPhivUn.ref(),
                    phaseSystemPtr_->phases()[phasei].phiv(),
                    alpha
                );
            }
            else
            {
                correctInflowOutflow(alphaPhivUn.ref(), phiv, alpha);
            }

            autoPtr<volScalarField::Internal> Su, Sp;
            calcSuSp(phasei, Su, Sp);

            if (MULESCorr_)
            {
                alphaPhivCorrs.set
                (
                    phasei,
                    alphaPhivUn() - alphaPhivs_[phasei]
                );
                alphaPhivUns.set(phasei, alphaPhivUn);
                surfaceScalarField& alphaPhivCorr = alphaPhivCorrs[phasei];

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);
                    MULES::limitCorr
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhivUns[phasei],
                        alphaPhivCorr,
                        Sp(),
                        (-Sp()*alpha)(),
                        1,
                        0,
                        MULEScontrols
                    );
                }
                else
                {
                    MULES::limitCorr
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhivUns[phasei],
                        alphaPhivCorr,
                        Sp(),
                        (-Sp()*alpha)(),
                        1,
                        0,
                        MULEScontrols
                    );
                }
            }
            else
            {
                alphaPhivCorrs.set(phasei, alphaPhivUn.ptr());
                surfaceScalarField& alphaPhivCorr = alphaPhivCorrs[phasei];

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    // NOTE: This converts alphaPhivCorr from a flux to a correction flux
                    MULES::limit
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        phiv,
                        alphaPhivCorr,
                        Sp(),
                        (Su() + divU*min(alpha.internalField(), scalar(1)))(),
                        1,
                        0,
                        true,
                        MULEScontrols
                    );
                }
                else
                {
                    // NOTE: This converts alphaPhivCorr from a flux to a correction flux
                    MULES::limit
                    (
                        1.0/mesh_.time().deltaT().value(),
                        geometricOneField(),
                        alpha,
                        phiv,
                        alphaPhivCorr,
                        Sp(),
                        (Su() + divU*min(alpha.internalField(), scalar(1)))(),
                        1,
                        0,
                        true,
                        MULEScontrols
                    );
                }
            }
        }

        if (passiveIndex_ < 0 || phases.size() > 2)
        {
            // We still need to limit the sum even when there is a passive phase
            // if there are more than 2 phases, as there could be interfaces
            // where the passive phase is absent and can't absorb a discrepancy

            MULES::limitSum(alphaPhivCorrs);
        }

        volScalarField sumAlpha
        (
            "sumAlpha",
            phiv.db(),
            mesh_,
            dimless,
            0
        );

        forAll(phases, phasei)
        {
            if (phasei == passiveIndex_)
            {
                continue;
            }

            volScalarField& alpha = thermoPtr_->alphas()[phasei];

            autoPtr<volScalarField::Internal> Su, Sp;
            calcSuSp(phasei, Su, Sp);

            if (MULESCorr_)
            {
                tmp<volScalarField::Internal> Spu = -Sp()*alpha;
                if (aCorr)
                {
                    Spu.ref() += (divU*(alpha.internalField()-alpha0[phasei].internalField()))();
                }
                if (nAlphaCorr_)
                {
                    alpha0.set(phasei, new volScalarField("alpha0", alpha));
                }

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    MULES::correct
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhivUns[phasei],
                        alphaPhivCorrs[phasei],
                        Sp(),
                        Spu()
                    );
                }
                else
                {
                    MULES::correct
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhivUns[phasei],
                        alphaPhivCorrs[phasei],
                        Sp(),
                        Spu()
                    );
                }

                // Under-relax the correction for all but the 1st corrector
                if (aCorr == 0)
                {
                    alphaPhivs_[phasei] += alphaPhivCorrs[phasei];
                }
                else
                {
                    alpha = 0.5*alpha + 0.5*alpha0[phasei];
                    alpha.correctBoundaryConditions();
                    alphaPhivs_[phasei] += 0.5*alphaPhivCorrs[phasei];
                }

            }
            else
            {
                tmp<surfaceScalarField> talphaPhiv =
                    upwind<scalar>(mesh_, phiv).flux(alpha)
                  + alphaPhivCorrs[phasei];

                if (phaseSystemPtr_)
                {
                    correctInflowOutflow
                    (
                        talphaPhiv.ref(),
                        phaseSystemPtr_->phases()[phasei].phiv(),
                        alpha
                    );
                }

                alphaPhivs_[phasei] = talphaPhiv;

                if (fv::localEulerDdt::enabled(mesh_))
                {
                    const volScalarField& rDeltaT =
                        fv::localEulerDdt::localRDeltaT(mesh_);

                    MULES::explicitSolve
                    (
                        rDeltaT,
                        geometricOneField(),
                        alpha,
                        alphaPhivs_[phasei],
                        Sp(),
                        // Divergence term is handled explicitly to be
                        // consistent with the explicit transport solution
                        (Su() + divU*min(alpha.internalField(), scalar(1)))()
                    );
                }
                else
                {
                    MULES::explicitSolve
                    (
                        1.0/mesh_.time().deltaTValue(),
                        geometricOneField(),
                        alpha,
                        alphaPhivs_[phasei],
                        Sp(),
                        // Divergence term is handled explicitly to be
                        // consistent with the explicit transport solution
                        (Su() + divU*min(alpha.internalField(), scalar(1)))()
                    );
                }
            }

            if (passiveIndex_ != -1 || aCorr == nAlphaCorr_-1)
            {
                sumAlpha += alpha;
            }

            if (aCorr == nAlphaCorr_-1)
            {
                Info<< alpha.group() << " volume fraction, min, max = "
                    << alpha.weightedAverage(mesh_.V()).value()
                    << ' ' << min(alpha).value()
                    << ' ' << max(alpha).value()
                    << endl;
            }
        }

        if (passiveIndex_ != -1)
        {
            volScalarField& alpha = thermoPtr_->alphas()[passiveIndex_];
            alpha.forceAssign(1-sumAlpha);
        }
        else if (aCorr == nAlphaCorr_-1)
        {
            Info<< "Phase-sum volume fraction, min, max = "
                << sumAlpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(sumAlpha).value()
                << ' ' << max(sumAlpha).value()
                << endl;
        }
    }

    autoPtr<surfaceScalarField> newPhi;
    if (phiPtr_)
    {
        newPhi.set
        (
            new surfaceScalarField
            (
                "newPhi",
                phiPtr_->db(),
                mesh_,
                phiPtr_->dimensions(),
                0
            )
        );
    }
    autoPtr<surfaceScalarField> sumAlphaPhiv;
    if (passiveIndex_ != -1)
    {
        sumAlphaPhiv.set
        (
            new surfaceScalarField
            (
                "sumAlphaPhiv",
                phiv.db(),
                mesh_,
                phiv.dimensions(),
                0
            )
        );
    }

    forAll(phases, phasei)
    {
        if (phasei == passiveIndex_)
        {
            continue;
        }

        volScalarField& alpha = thermoPtr_->alphas()[phasei];

        surfaceScalarField& alphaPhiv = alphaPhivs_[phasei];

        if (alphaApplyPrevCorr_)
        {
            alphaPhivCorr0_.set
            (
                phasei,
                alphaPhiv - alphaPhivCorr0_[phasei]
            );
            alphaPhivCorr0_[phasei].rename(alpha.name()+"PhivCorr0");
        }

        if (crankNicolson_)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhiv =
                (alphaPhiv - (1.0 - cnCoeff_)*alphaPhiv.oldTime())/cnCoeff_;
        }

        if (phiPtr_)
        {
            newPhi() +=
                fvc::interpolate(thermoPtr_->thermos()[phasei].rho())*alphaPhiv;
        }

        if (passiveIndex_ != -1)
        {
            sumAlphaPhiv() += alphaPhiv;
        }
    }

    if (passiveIndex_ != -1)
    {
        alphaPhivs_[passiveIndex_].forceAssign(phiv - sumAlphaPhiv());

        if (phiPtr_)
        {
            // Calculate the end-of-time-step mass flux (phiv rather than phiCN)
            newPhi() +=
                fvc::interpolate(thermoPtr_->thermos()[passiveIndex_].rho())
               *(phiv-sumAlphaPhiv());
        }
    }

    if (phiPtr_)
    {
        *phiPtr_ = newPhi();
    }

    thermoPtr_->fractions().recomputeCombinedAlphas();
}


void Foam::fv::MULESVolumeFractionSolver::correctInflowOutflow
(
    surfaceScalarField& alphaPhi,
    const surfaceScalarField& phiv,
    const volScalarField& alpha
) const
{
    if (correctInflowOutflow_)
    {
        surfaceScalarField::Boundary& alphaPhivBf = alphaPhi.boundaryFieldRef();
        const volScalarField::Boundary& alphaBf = alpha.boundaryField();

        const surfaceScalarField::Boundary& phivBf = phiv.boundaryField();

        forAll(alphaPhivBf, patchi)
        {
            fvsPatchScalarField& alphaPhivp = alphaPhivBf[patchi];

            if (!alphaPhivp.coupled())
            {
                alphaPhivp = phivBf[patchi]*alphaBf[patchi];
            }
        }
    }
}


void Foam::fv::MULESVolumeFractionSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    const Time& runTime = mesh().time();

    PtrList<volScalarField>& alphas = thermoPtr_->alphas();

    volScalarField::Internal divU
    (
        phaseSystemPtr_
        // U here is only used to select the correct time scheme
      ? fvc::div(fvc::absolute(*phivPtr_, phaseSystemPtr_->phases()[0].U()))
      : fvc::div(fvc::absolute(*phivPtr_, *UPtr_))
    );

    if (nAlphaSubCycles_ > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        tmp<volScalarField> trSubDeltaT;
        if (LTS_)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles_);
        }

        autoPtr<surfaceScalarField> phiSum;
        PtrList<surfaceScalarField> alphaPhivSums(alphas.size());
        if (phiPtr_)
        {
            phiSum.set
            (
                new surfaceScalarField
                (
                    "phiSum",
                    alphas[0].db(),
                    mesh_,
                    dimMass/dimTime,
                    0
                )
            );
        }
        forAll(alphaPhivSums, phasei)
        {
            alphaPhivSums.set
            (
                phasei,
                new surfaceScalarField
                (
                    "phivSum" + alphas[phasei].name(),
                    alphas[0].db(),
                    mesh_,
                    dimVolume/dimTime,
                    0
                )
            );
        }

        // Stores old-times for later restoration. This must come before
        // sub-cycling is started through creation of subCycleTime below
        PtrList<subCycleField<volScalarField>> alphaSubCycles(alphas.size());
        forAll(alphas, phasei)
        {
            alphaSubCycles.set
            (
                phasei,
                new subCycleField<volScalarField>(alphas[phasei])
            );
        }

        autoPtr<subCycleTime> subCycle
        (
            new subCycleTime
            (
                const_cast<Time&>(mesh_.time()),
                nAlphaSubCycles_
            )
        );

        forAll(alphas, phasei)
        {
            alphaSubCycles[phasei].updateTimeIndex();
        }

        while (!(++subCycle()).end())
        {
            solveAlphas(divU);
            if (phiPtr_)
            {
                phiSum() += (runTime.deltaT()/totalDeltaT)*(*phiPtr_);
            }
            forAll(alphas, phasei)
            {
                alphaPhivSums[phasei] += alphaPhivs_[phasei];
            }
        }

        // End sub-cycling
        subCycle.clear();
        // Restore old-time fields
        alphaSubCycles.clear();

        if (phiPtr_)
        {
            *phiPtr_ = phiSum;
        }
        forAll(alphas, phasei)
        {
            alphaPhivs_[phasei] = alphaPhivSums[phasei]/nAlphaSubCycles_;
        }
    }
    else
    {
        solveAlphas(divU);
    }

    thermoPtr_->correct();
}


void Foam::fv::MULESVolumeFractionSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName && finalIter)
    {
        Info<< endl;

        if (!LTS_)
        {
            // Calculate Courant number at interface

            maxAlphaCo_ = 0.0;
            scalar meanAlphaCoNum = 0.0;

            if (mesh().nInternalFaces())
            {
                scalarField sumPhi
                (
                    thermoPtr_->nearInterface()().primitiveField()
                   *fvc::surfaceSum(mag(*phivPtr_))().primitiveField()
                );

                maxAlphaCo_ =
                    0.5*gMax(sumPhi/mesh().V().field())*mesh().time().deltaTValue();

                meanAlphaCoNum =
                    0.5
                   *(
                        gSum(sumPhi)/gSum(mesh().V().field())
                    )*mesh().time().deltaTValue();
            }

            Info<< "Interface Courant Number mean: " << meanAlphaCoNum
                << " max: " << maxAlphaCo_ << endl;
        }
    }
}


void Foam::fv::MULESVolumeFractionSolver::topoChange(const polyTopoChangeMap& map)
{
    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    alphaPhivCorr0_.clear();
}


bool Foam::fv::MULESVolumeFractionSolver::isDeprecatedAlphaControl
(
    const dictionary& alphaControls, const word& key
)
{
    if (alphaControls.found(key))
    {
        DeprecationIOWarningInFunction
        (
            alphaControls,
            key,
            "keyword",
            40300,
            "Please specify alpha controls in the fvOptions."
          + this->name() + " subdictionary instead of in fvSolution"
        );
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::fv::MULESVolumeFractionSolver::setRDeltaT()
{
    if (LTS_)
    {
        const dictionary& solnDict = solnControlPtr_->dict();

        const scalar maxCo = solnDict.lookupOrDefault<scalar>("maxCo", 0.9);
        const scalar maxAlphaCo =
            solnDict.lookupOrDefault<scalar>("maxAlphaCo", 0.2);
        const scalar rDeltaTSmoothingCoeff =
            solnDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1);
        const label nAlphaSpreadIter =
            solnDict.lookupOrDefault<label>("nAlphaSpreadIter", 1);
        const scalar alphaSpreadDiff =
            solnDict.lookupOrDefault<scalar>("alphaSpreadDiff", 0.2);
        const scalar alphaSpreadMax =
            solnDict.lookupOrDefault<scalar>("alphaSpreadMax", 0.99);
        const scalar alphaSpreadMin =
            solnDict.lookupOrDefault<scalar>("alphaSpreadMin", 0.01);
        const label nAlphaSweepIter =
            solnDict.lookupOrDefault<label>("nAlphaSweepIter", 5);
        const scalar rDeltaTDampingCoeff =
            solnDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0);
        const scalar maxDeltaT =
            solnDict.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        volScalarField& rDeltaT = lookupOrConstructRDeltaT();

        // Set the reciprocal time-step from the local Courant number
        tmp<volScalarField::Internal> rDeltaTT
        (
            volScalarField::Internal::New
            (
                typeName + "-rDeltaT",
                fvc::surfaceSum(mag(*phivPtr_))()()/((2*maxCo)*mesh_.V())
            )
        );

        rDeltaT.ref() = max(rDeltaT, rDeltaTT());

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        // Create a compound 'interface' indicater that it the product
        // of alphas
        volScalarField alpha(thermoPtr_->alphas()[0]);
        for
        (
            label phasei = 1;
            phasei < thermoPtr_->alphas().size();
            phasei++
        )
        {
            alpha *= thermoPtr_->alphas()[phasei];
        }

        if (maxAlphaCo < maxCo)
        {
            // Further limit the reciprocal time-step
            // in the vicinity of the interface
            volScalarField alphaBar(fvc::average(alpha));

            rDeltaT.ref() =
                max
                (
                    rDeltaT,
                    pos0(alphaBar() - alphaSpreadMin)
                   *pos0(alphaSpreadMax - alphaBar())
                   *fvc::surfaceSum(mag(*phivPtr_))()()
                   /((2*maxAlphaCo)*mesh_.V())
                );
        }

        // Update the boundary values of the reciprocal time-step
        rDeltaT.correctBoundaryConditions();

        if (rDeltaTSmoothingCoeff < 1.0)
        {
            fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
        }

        if (nAlphaSpreadIter > 0)
        {
            fvc::spread
            (
                rDeltaT,
                alpha,
                nAlphaSpreadIter,
                alphaSpreadDiff,
                alphaSpreadMax,
                alphaSpreadMin
            );
        }

        if (nAlphaSweepIter > 0)
        {
            fvc::sweep(rDeltaT, alpha, nAlphaSweepIter, alphaSpreadDiff);
        }

        // Limit rate of change of time scale
        // - reduce as much as required
        // - only increase at a fraction of old time scale
        if
        (
            rDeltaTDampingCoeff < 1.0
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

        Info<< "    Flow = "
            << gMin(1/rDeltaT.primitiveField()) << ", "
            << gMax(1/rDeltaT.primitiveField()) << endl;

        return true;
    }
    return false;
}

// ************************************************************************* //
