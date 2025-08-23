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
    (c) 2019-2025 Engys Ltd.
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fluidEnergySolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "interpolation/surfaceInterpolation/schemes/midPoint/midPoint.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "finiteVolume/ddtSchemes/steadyStateDdtScheme/steadyStateDdtScheme.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "equationOfState/BoussinesqLaw/BoussinesqLaw.H"
#include "equationOfState/rhoConstLaw/rhoConstLaw.H"
#include "basicThermo/BasicThermo.H"
#include "eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "eulerianMultiphaseSystem/eulerianMultiphaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fluidEnergySolver, 0);
}
}

makeFvSolverOption(fluidEnergySolver);

const Foam::Enum<Foam::fv::fluidEnergySolver::formulationType>
Foam::fv::fluidEnergySolver::formulationTypeNames_
(
    {
        {formulationType::totalEnergyEnergy, "totalEnergyEnergy"},
        {formulationType::totalEnergyTemperature, "totalEnergyTemperature"},
        {formulationType::boussinesqEnergy, "boussinesqEnergy"},
        {
            formulationType::boussinesqTemperature,
            "boussinesqTemperature"
        },
        {formulationType::phasicEnergy, "phasicEnergy"}
    }
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void Foam::fv::fluidEnergySolver::createK()
{
    if (!K_.valid())
    {
        word UName = addPhaseName("U");
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(UName);
        K_.set(new volScalarField("K", 0.5*magSqr(U)));

        if (U.nOldTimes())
        {
            const volVectorField* Uold = &U.oldTime();
            volScalarField* Kold = &K_().oldTime();
            Kold->forceAssign(0.5*magSqr(*Uold));

            while (Uold->nOldTimes())
            {
                Uold = &Uold->oldTime();
                Kold = &Kold->oldTime();
                Kold->forceAssign(0.5*magSqr(*Uold));
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::fv::fluidEnergySolver::divMRF
(
    const surfaceScalarField& phi,
    const volScalarField* alpha
)
{
    const volScalarField& p = thermoPtr_->p();
    tmp<surfaceScalarField> tphiMRF
    (
        surfaceScalarField::New
        (
            "phiMRF",
            mesh_,
            dimensionedScalar(phi.dimensions()/dimDensity, 0),
            calculatedFvsPatchField<scalar>::typeName
        )
    );
    fvOptions().makeAbsolute(tphiMRF.ref());

    if (alpha)
    {
        tphiMRF.ref() *= fvc::interpolate(*alpha);
    }

    return fv::gaussConvectionScheme<scalar>
    (
        mesh_,
        tphiMRF(),
        linear<scalar>(mesh_)
    ).fvcDiv(tphiMRF(), p);
}


tmp<Foam::volScalarField> Foam::fv::fluidEnergySolver::viscousStressTerm
(
    const volVectorField& U
) const
{
    const compressible::turbulenceModel& turbulence =
        obr_.lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                phaseName_
            )
        );
    volScalarField muEff("muEff", turbulence.muEff());
    volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));
    surfaceScalarField sigmaDotU
    (
        "sigmaDotU",
        (
            fvc::interpolate(muEff)*mesh_.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh_.Sf(), tauMC)
        )
      & (linearInterpolate(U))
    );
    return fvc::div(sigmaDotU);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fluidEnergySolver::fluidEnergySolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    thermoPtr_(nullptr),
    phaseSystemPtr_(nullptr),
    hasMRF_(false),
    lowFroudeStabilisation_
    (
        dict.lookupOrDefault<Switch>("lowFroudeStabilisation", true)
    ),
    lowFroudeStabilisationCoeff_
    (
        dict.lookupOrDefault<scalar>("lowFroudeStabilisationCoeff", 0.05)
    ),
    viscousStressTerm_
    (
        dict.lookupOrDefault<Switch>("viscousStressTerm", false)
    ),
    formulation_(totalEnergyEnergy)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::fv::fluidEnergySolver::fvmDiv
(
    const surfaceScalarField& phi,
    const volScalarField& he
)
{
    // Use multivarible convection scheme if available
    typedef speciesConcentrationSolver speciesType;
    const wordList speciesTransportNames(obr_.names<speciesType>());
    if (speciesTransportNames.size() == 1)
    {
        const word& speciesTransportName = speciesTransportNames.first();
        speciesType& concentrationTransport =
            obr_.lookupObjectRef<speciesType>(speciesTransportName);
        return concentrationTransport.convTerm().fvmDiv(phi, he);
    }

    return fvm::div(phi, he);
}


void Foam::fv::fluidEnergySolver::addEnergyToMultivariableScheme()
{
    const wordList speciesTransportNames(obr_.names<speciesType>());
    if
    (
        speciesTransportNames.size() == 1
     && (formulation_ == boussinesqEnergy || formulation_ == totalEnergyEnergy)
    )
    {
        const word& speciesTransportName = speciesTransportNames.first();
        obr_.lookupObjectRef<speciesType>(speciesTransportName).addField
        (
            thermoPtr_->he()
        );
    }
}


tmp<volScalarField> Foam::fv::fluidEnergySolver::lookupOrConstructDpdt
(
    const volScalarField& p
)
{
    //- Grab ddt(p) based on energy ddt.
    //  If momentum-continuity and energy is setup in different scales,
    //  ddt(p) in continuity and ddt(p) here must be different.
    //  Therefore this one is based on ddt(rho, h/T)
    tmp<ddtScheme<scalar>> ddtpscheme = fv::getDdtScheme(p);
    {
        basicThermo& thermo = *thermoPtr_;
        if
        (
            (formulation_ == totalEnergyEnergy)
          || (formulation_ == boussinesqEnergy)
        )
        {
            const volScalarField& he = thermo.he();
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(addPhaseName("rho"));
            ddtpscheme = fv::getDdtScheme(rho,he);
        }
        else if
        (
            (formulation_ == totalEnergyTemperature)
          || (formulation_ == boussinesqTemperature)
        )
        {
            const volScalarField& T = thermo.T();
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(addPhaseName("rho"));
            ddtpscheme = fv::getDdtScheme(rho,T);
        }
    }

    // TODO: Could be removed after testing of new solution order
    return
        isReactingFoam()
      ? volScalarField::New
        (
            "dpdt",
            obr_.lookupObject<volScalarField>("dpdt")
        )
      : ddtpscheme->fvcDdt(p);
}


bool Foam::fv::fluidEnergySolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    thermoPtr_ = &multiphaseThermo::lookupOrCreate(obr(), phaseName_);

    if (dict_.found("formulation"))
    {
        formulation_ = formulationTypeNames_.lookup("formulation", dict_);
    }
    else
    {
        // Check for eulerianPhaseSystem
        if (obr_.foundObject<eulerianPhaseSystem>(eulerianPhaseSystem::typeName))
        {
            formulation_ = phasicEnergy;
        }
        else if (isA<multiphaseThermo>(*thermoPtr_))
        {
            formulation_ = totalEnergyTemperature;
        }
        else if
        (
            thermoPtr_->distinctBuoyancy()
         && thermoPtr_->he().member() == "e"
         && thermoPtr_->materials().sTable
            (
                word::null, word::null
            )["buoyancy"]->type() == BoussinesqLaw::typeName
         && thermoPtr_->materials().sTable
            (
                word::null, word::null
            )[rhoModel::typeName]->type() == rhoConstLaw::typeName
        )
        {
            formulation_ = boussinesqEnergy;
        }
        else
        {
            formulation_ = totalEnergyEnergy;
        }
    }

    if (formulation_ == totalEnergyEnergy)
    {
        Info<< "Solving for internal energy in total energy formulation"
            << endl;
        thermoPtr_->setToCalculateTFromhe(true);
    }
    else if (formulation_ == totalEnergyTemperature)
    {
        Info<< "Solving for temperature in total energy formulation" << endl;
        thermoPtr_->setToCalculateTFromhe(false);
    }
    else if (formulation_ == boussinesqEnergy)
    {
        Info<< "Solving for internal energy in Boussinesq formulation"
            << endl;
        thermoPtr_->setToCalculateTFromhe(true);
    }
    else if (formulation_ == boussinesqTemperature)
    {
        Info<< "Solving for temperature in Boussinesq formulation" << endl;
        thermoPtr_->setToCalculateTFromhe(false);
    }
    else if (formulation_ == phasicEnergy)
    {
        Info<< "Solving for energy in phase system" << endl;
        phaseSystemPtr_ =
            &refCast<eulerianMultiphaseSystem>
            (
                obr_.lookupObjectRef<eulerianPhaseSystem>
                (
                    eulerianPhaseSystem::typeName
                )
            );
        thermoPtr_->setToCalculateTFromhe(true);
    }
    else
    {
        NotImplemented;
    }

    if (formulation_ == phasicEnergy && !thermoPtr_->isPhasicVariable("T"))
    {
        FatalErrorInFunction
            << "Phase-specific T fields are required for this solver, "
            << "but were not found."
            << exit(FatalError);
    }
    else if (formulation_ != phasicEnergy && thermoPtr_->isPhasicVariable("T"))
    {
        FatalErrorInFunction
            << "Phase-specific T fields were found but are not used by this "
            << "solver."
            << exit(FatalError);
    }

    addEnergyToMultivariableScheme();

    mesh_.schemes().setFluxRequired(thermoPtr_->he().name());
    mesh_.schemes().setFluxRequired(thermoPtr_->T().name());

    // Check if we have MRF
    hasMRF_ = false;
    for (int i=0; i < fvOptions().optionList::size(); ++i)
    {
        if (fvOptions().optionList::operator[](i).isMRF())
        {
            hasMRF_ = true;
        }
    }

    if
    (
        formulation_ == totalEnergyEnergy
     || formulation_ == totalEnergyTemperature
    )
    {
        Info<< "Creating field K\n" << endl;
        createK();
    }

    return true;
}


void Foam::fv::fluidEnergySolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    if (formulation_ == phasicEnergy)
    {
        eulerianMultiphaseSystem::phaseModelList& phases =
            phaseSystemPtr_->phases();

        DynamicList<word> UNames, heNames;
        forAll(phases, phasei)
        {
            eulerianPhaseModel& phase = phases[phasei];
            UNames.append(phase.U()().name());
            heNames.append(phase.thermo().he().name());
        }

        forAll(phases, phasei)
        {
            solveNames.append(heNames[phasei]);
            requiredDependencies.insert(heNames[phasei], UNames);
            optionalDependencies.insert(heNames[phasei], {"fvMesh"});
        }

        // Energy correctors
        correctorMembers.insert("energyCorrector", heNames);
    }
    else
    {
        // Solves
        word heTName = thermoPtr_->heT().name();
        solveNames.append({heTName});

        // Dependencies
        if (!isReactingFoam())
        {
            optionalDependencies.insert
            (
                heTName,
                {addPhaseName("U"), thermoPtr_->p().name(), "fvMesh"}
            );
        }
    }

    optionalDependencies.insert("turbulence", solveNames);

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


bool Foam::fv::fluidEnergySolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        nEnergyCorrectors_ =
            solnControlPtr_->dict().lookupOrDefault<label>
            (
                "nEnergyCorrectors", 1
            );
        return false;
    }
    else if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if (correctorName == "energyCorrector")
    {
        // Update nut/alphat
        phaseSystemPtr_->correctEnergyTransport();
        return (corrector >= nEnergyCorrectors_-1);
    }
    else
    {
        return false;
    }
}


tmp<fvScalarMatrix>
Foam::fv::fluidEnergySolver::assembleScalarMatrix
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    basicThermo& thermo = *thermoPtr_;
    if (formulation_ == totalEnergyEnergy)
    {
        // Conservative formulation written in terms of energy
        // Typically used in single-phase applications

        volScalarField& he = thermo.he();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));

        // Update turbulence kinetic energy
        K_() = 0.5*magSqr(U);
        const volScalarField& K = K_();

        const volScalarField& p = thermo.p();
        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, he)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn
        (
            fvmDiv(phi, he) + fvc::div(phi, K)
         ==
            fvOptions()(rho, he)
        );

        if (viscousStressTerm_)
        {
            EEqn.ref() -= this->viscousStressTerm(U);
        }

        // Pick up any fvOptions in T to provide a consistent interface where
        // sources can be specified in terms of T regardless of formulation
        volScalarField& T = thermo.T();
        if (!isReactingFoam())
        {
            EEqn.ref() -=
                changeVariable
                (
                    fvOptions()(rho*thermo.Cv(), T),
                    1/thermo.Cpv(),
                    he
                );
        }
        // updateCoeffs was called in matrix construction but evaluate was
        // not called - resetUpdate reverts to updatable state, to avoid
        // interfering with normal update of boundaries
        forAll(T.boundaryField(), patchi)
        {
            T.boundaryFieldRef()[patchi].resetUpdate();
        }

        if (!steady)
        {
            EEqn.ref() += fvm::ddt(rho, he) + fvc::ddt(rho, K);
        }

        if (he.member() == "e")
        {
            EEqn.ref() +=
                fvc::div
                (
                    fvc::absolute(phi, rho, U), thermo.pAbs()/rho, "div(phiv,p)"
                );
        }
        else if (!steady && thermoPtr_->dpdt())
        {
            // Update pressure time derivative if needed
            volScalarField dpdt(lookupOrConstructDpdt(p));
            if (mesh_.moving())
            {
                dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
            }
            EEqn.ref() -= dpdt;
        }

        // Add MRF contribution
        if (hasMRF_)
        {
            EEqn.ref() += divMRF(phi);
        }

        // Add the term laplacian(kappaEff, T) to RHS
        if (thermo.isCpvConst() || isReactingFoam())
        {
            // Faster calculation
            tmp<volScalarField> alphaEff(turbulence.alphaEff());
            word scheme = "laplacian(kappa";
            if (alphaEff().group() != word::null)
            {
                scheme += "." + alphaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -= fvm::laplacian(alphaEff(), he, scheme);
        }
        else
        {
            tmp<volScalarField> kappaEff(turbulence.kappaEff());
            word scheme = "laplacian(kappa";
            if (kappaEff().group() != word::null)
            {
                scheme += "." + kappaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -=
                changeVariable
                (
                    fvm::laplacian(kappaEff(), thermo.T(), scheme),
                    1/thermo.Cpv(),
                    he
                );
        }

        if (thermo.buoyant())
        {
            const uniformDimensionedVectorField& g =
                obr_.lookupObject<uniformDimensionedVectorField>("g");
            tmp<volScalarField> bRho
            (
                thermo.distinctBuoyancy()
              ? thermo.buoyantRho()
              : tmp<volScalarField>(rho)
            );
            EEqn.ref() -= bRho()*(U & g);

            if (lowFroudeStabilisation_ && lowFroudeStabilisationCoeff_ > 0)
            {
                // Add stabilisation to the diagonal for the above term
                dimensionedScalar gCmptSum
                (
                    "gCmptSum", g.dimensions(), cmptSum(cmptMag(g.value()))
                );
                tmp<volScalarField> offDiagCoeff
                (
                    lowFroudeStabilisationCoeff_*bRho*gCmptSum
                );
                offDiagCoeff->dimensions() *= U.dimensions()/he.dimensions();

                if (debug)
                {
                    volScalarField D("D", U.db(), mesh_, dimless);
                    D.primitiveFieldRef() = EEqn().D()()/mesh_.V();
                    volScalarField stabCoeff("stabCoeff", U.db(), mesh_, dimless);
                    stabCoeff.primitiveFieldRef() = offDiagCoeff();
                    if (mesh_.time().writeTime())
                    {
                        D.write();
                        stabCoeff.write();
                    }
                }

                // offDiagCoeff is the diagonal dominance that would be required
                // for these coupled terms. Ensure diagonal is at least as large
                offDiagCoeff->primitiveFieldRef() -=
                    min
                    (
                        max(EEqn().D()()/mesh_.V(), scalar(0)),
                        offDiagCoeff().primitiveField()
                    );
                EEqn.ref() += fvm::SuSp(offDiagCoeff(), he);
                EEqn.ref() -= fvc::Sp(offDiagCoeff, he);
            }
        }

        if (thermoPtr_->properties().found("porosity"))
        {
            EEqn.ref() *= thermoPtr_->properties().lookup<scalar>("porosity");
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == totalEnergyTemperature)
    {
        // Conservative formulation written in terms of temperature. Typically
        // used for multiphase to avoid large discontinuities in energy.

        volScalarField& T = thermo.T();
        volScalarField& rho =
            obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const surfaceScalarField& phiv =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phiv"));

        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));

        // Update turbulence kinetic energy
        K_() = 0.5*magSqr(U);
        const volScalarField& K = K_();

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, T)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<volScalarField> Cv;
        if (isA<multiphaseThermo>(*thermoPtr_))
        {
            // Perform harmonic averaging of Cv consistent with averaging
            // temperature rather than energy
            Cv = 1/refCast<multiphaseThermo>(*thermoPtr_).rCv();
        }
        else
        {
            Cv = volScalarField::New(thermoPtr_->Cv().name(), thermoPtr_->Cv());
        }
        tmp<fvScalarMatrix> EEqn
        (
            Cv()*(fvm::div(phi, T) + (thermo.TRef()*fvc::div(phi)))
          + (
                fvc::div
                (
                    fvc::absolute(phiv, U), thermo.pAbs(), "div(phiv,p)"
                )
              + fvc::ddt(rho, K) + fvc::div(phi, K)
            )
         == fvOptions()(rho*Cv(), T)
        );
        if (!steady)
        {
            EEqn.ref() += Cv*(fvm::ddt(rho, T) + thermo.TRef()*fvc::ddt(rho));
        }

        // Add MRF contribution
        if (hasMRF_)
        {
            EEqn.ref() += divMRF(phi);
        }

        // Add the term laplacian(kappaEff, T) to RHS
        tmp<volScalarField> kappaEff(turbulence.kappaEff());
        word scheme = "laplacian(kappa";
        if (kappaEff().group() != word::null)
        {
            scheme += "." + kappaEff().group();
        }
        scheme += "," + thermo.T().name() + ")";
        EEqn.ref() -=
            fvm::laplacian(kappaEff(), thermo.T(), scheme);

        if (thermoPtr_->properties().found("porosity"))
        {
            EEqn.ref() *= thermoPtr_->properties().lookup<scalar>("porosity");
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == boussinesqEnergy)
    {
        // Not conservative, and viscous dissipation terms omitted which are
        // typically negligible in Boussinesq approximation

        volScalarField& he = thermo.he();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));
        const volVectorField& U =
            obr_.lookupObject<volVectorField>(addPhaseName("U"));
        const volScalarField& p = thermo.p();

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, he)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn(fvmDiv(phi, he));

        if (!steady)
        {
            EEqn.ref() += fvm::ddt(rho, he);
        }

        if (he.name() == "h")
        {
            EEqn.ref() -=
                fvc::div
                (
                    fvc::absolute(phi, rho, U),
                    thermo.pAbs()/rho,
                    "div(phiv,p)"
                );
            if (!steady && thermoPtr_->dpdt())
            {
                volScalarField dpdt(fvc::ddt(p));
                if (mesh_.moving())
                {
                    dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
                }
                EEqn.ref() -= dpdt;
            }

            // Add MRF contribution
            if (hasMRF_)
            {
                EEqn.ref() += divMRF(phi);
            }
        }

        EEqn.ref() -= fvOptions()(rho, he);

        // Pick up any fvOptions in T to provide a consistent interface where
        // sources can be specified in terms of T regardless of formulation
        volScalarField& T = thermo.T();
        EEqn.ref() -=
            changeVariable
            (
                fvOptions()(rho*thermo.Cp(), T),
                1/thermo.Cpv(),
                he
            );
        // updateCoeffs was called in matrix construction but evaluate was
        // not called - resetUpdate reverts to updatable state, to avoid
        // interfering with normal update of boundaries
        forAll(T.boundaryField(), patchi)
        {
            T.boundaryFieldRef()[patchi].resetUpdate();
        }


        // Add the term laplacian(kappaEff, T) to RHS
        if (thermo.isCpvConst())
        {
            // Faster calculation
            tmp<volScalarField> alphaEff(turbulence.alphaEff());
            word scheme = "laplacian(kappa";
            if (alphaEff().group() != word::null)
            {
                scheme += "." + alphaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -= fvm::laplacian(alphaEff(), he, scheme);
        }
        else
        {
            tmp<volScalarField> kappaEff(turbulence.kappaEff());
            word scheme = "laplacian(kappa";
            if (kappaEff().group() != word::null)
            {
                scheme += "." + kappaEff().group();
            }
            scheme += "," + thermo.T().name() + ")";
            EEqn.ref() -=
                changeVariable
                (
                    fvm::laplacian(kappaEff(), thermo.T(), scheme),
                    1/thermo.Cpv(),
                    he
                );
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == boussinesqTemperature)
    {
        // Not conservative, and visous dissipation terms omitted which are
        // typically negligible in Boussinesq approximation

        volScalarField& T = thermo.T();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(addPhaseName("phi"));

        bool steady =
            isA<steadyStateDdtScheme<scalar>>(fv::getDdtScheme(rho, T)());

        const compressible::turbulenceModel& turbulence =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                IOobject::groupName
                (
                    compressible::turbulenceModel::propertiesName,
                    phaseName_
                )
            );

        tmp<fvScalarMatrix> EEqn
        (
            fvm::div(phi, T) + (thermo.TRef()*fvc::div(phi))
        );

        if (!steady)
        {
            EEqn.ref() += (fvm::ddt(rho, T) + thermo.TRef()*fvc::ddt(rho));
        }

        tmp<volScalarField> Cv = thermo.Cv();
        EEqn.ref() *= Cv();
        EEqn.ref() -= fvOptions()(rho*Cv, T);

        // Add the term laplacian(kappaEff, T) to RHS
        tmp<volScalarField> kappaEff(turbulence.kappaEff());
        word scheme = "laplacian(kappa";
        if (kappaEff().group() != word::null)
        {
            scheme += "." + kappaEff().group();
        }
        scheme += "," + thermo.T().name() + ")";
        EEqn.ref() -= fvm::laplacian(kappaEff(), thermo.T(), scheme);

        if (thermoPtr_->properties().found("porosity"))
        {
            EEqn.ref() *= thermoPtr_->properties().lookup<scalar>("porosity");
        }

        EEqn->relax();

        fvOptions().constrain(EEqn.ref());

        return EEqn;
    }
    else if (formulation_ == phasicEnergy)
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
        if (solveMember == "h" || solveMember == "e")
        {
            eulerianPhaseModel& phase = phaseSystemPtr_->phases()[solveGroup];
            const label phasei = phase.index();
            return assemblePhasicEEqn(phasei);
        }
        else
        {
            return tmp<fvScalarMatrix>();
        }
    }
    else
    {
        NotImplemented;
    }
}


tmp<fvScalarMatrix>
Foam::fv::fluidEnergySolver::assemblePhasicEEqn
(
    const label phasei
)
{
    if (!heatTransferPtr_.valid())
    {
        heatTransferPtr_ = phaseSystemPtr_->heatTransfer();
    }
    eulerianPhaseSystem::heatTransferTable& heatTransfer = heatTransferPtr_();

    eulerianPhaseModel& phase = phaseSystemPtr_->phases()[phasei];

    const volScalarField& alpha = phase;
    tmp<volScalarField> trho(phase.rho());
    const volScalarField& rho = trho();
    const volVectorField& U = phase.U();

    tmp<fvScalarMatrix> EEqn(phase.heEqn());

    if (EEqn.valid())
    {
        // Add MRF contribution
        if (hasMRF_)
        {
            EEqn.ref() += divMRF(phase.alphaPhi(), &alpha);
        }

        EEqn.ref() -=
            *heatTransfer[phase.name()]
          + fvOptions()(alpha, rho, phase.thermo().he());

        // Pick up any fvOptions in T to provide a consistent interface where
        // sources can be specified in terms of T regardless of formulation
        volScalarField& T = phase.thermo().T();
        EEqn.ref() -=
            changeVariable
            (
                fvOptions()(rho*phase.thermo().Cp(), T),
                1/phase.thermo().Cpv(),
                phase.thermo().he()
            );
        // updateCoeffs was called in matrix construction but evaluate was
        // not called - resetUpdate reverts to updatable state, to avoid
        // interfering with normal update of boundaries
        forAll(T.boundaryField(), patchi)
        {
            T.boundaryFieldRef()[patchi].resetUpdate();
        }

        if (thermoPtr_->buoyant())
        {
            const uniformDimensionedVectorField& g =
                obr_.lookupObject<uniformDimensionedVectorField>("g");
            EEqn.ref() -= alpha*rho*(U & g);
        }

        EEqn->relax();
        fvOptions().constrain(EEqn.ref());
    }
    return EEqn;
}


void Foam::fv::fluidEnergySolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    if
    (
        formulation_ == totalEnergyEnergy
     || formulation_ == boussinesqEnergy
    )
    {
        // Solved for energy

        volScalarField& he = thermoPtr_->he();

        if (verbose_)
        {
            Info<< he.name() << " min/max: "
                << min(he).value() << " "
                << max(he).value() << endl;
        }

        fvOptions().correct(he);
    }
    else if
    (
        formulation_ == totalEnergyTemperature
     || formulation_ == boussinesqTemperature
    )
    {
        // Solved for temperature

        volScalarField& he = thermoPtr_->he();
        volScalarField& T = thermoPtr_->T();

        fvOptions().correct(T);

        he.forceAssign(thermoPtr_->he(thermoPtr_->p(), T));

        // In multiphase case, set energies of phase thermos such as to produce
        // the same temperature for all
        if (isA<multiphaseThermo>(*thermoPtr_))
        {
            multiphaseThermo& mpThermo = refCast<multiphaseThermo>(*thermoPtr_);
            forAll(mpThermo.thermos(), ti)
            {
                basicThermo& thermo = mpThermo.thermos()[ti];
                thermo.he().forceAssign(thermo.he(thermo.p(), T));
            }
        }
    }
    else if (formulation_ == phasicEnergy)
    {
    }
    else
    {
        NotImplemented;
    }

    if (!phaseSystemPtr_)
    {
        // Thermo correction is done in endIteration
        // if phaseSystemPtr_ is valid
        thermoPtr_->correct();
    }
}


void Foam::fv::fluidEnergySolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == "energyCorrector")
    {
        heatTransferPtr_.clear();

        phaseSystemPtr_->correctThermo();
        phaseSystemPtr_->correct();

        if (finalIter)
        {
            phaseSystemPtr_->correctKinematics();

            eulerianPhaseSystem::phaseModelList& phases =
                phaseSystemPtr_->phases();
            forAll(phases, phasei)
            {
                eulerianPhaseModel& phase = phases[phasei];

                Info<< phase.name() << " min/max T "
                    << min(phase.thermo().T()).value()
                    << " - "
                    << max(phase.thermo().T()).value()
                    << endl;
            }
        }
    }
    if (correctorName == solverObject::outerCorrectorName && finalIter)
    {
        // Update rho
        if (obr_.foundObject<volScalarField>(addPhaseName("rho")))
        {
            volScalarField& rho =
                obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));
            rho = thermoPtr_->rho();
            rho.relax();
        }
    }
}

// ************************************************************************* //
