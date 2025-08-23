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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "lagrangianSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "basicThermo/basicThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(lagrangianSolver, 0);
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::fv::lagrangianSolver::lagrangianTypes,
    7
>::names[] =
{
    "basicReactingCloud",
    "basicReactingMultiphaseCloud",
    "basicThermoCloud",
    "coalCloud",
    "basicSprayCloud",
    "basicKinematicCollidingCloud",
    "basicKinematicCloud"
};


const Foam::NamedEnum
<
    Foam::fv::lagrangianSolver::lagrangianSolver::lagrangianTypes,
    7
>
Foam::fv::lagrangianSolver::lagrangianSolver::lagrangianTypesNames_;


makeFvSolverOption(lagrangianSolver);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::SLGThermo& Foam::fv::lagrangianSolver::slgThermo()
{
    if (!this->mesh_.thisDb().found("SLGThermo"))
    {
        this->mesh_.thisDb().store(new SLGThermo(mesh_, *thermoPtr_));
    }
    return this->mesh_.thisDb().lookupObject<SLGThermo>("SLGThermo");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::lagrangianSolver::lagrangianSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{
    thermoPtr_ =
        &refCast<rhoThermo>(multiphaseThermo::lookupOrCreate(obr, phaseName_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::lagrangianSolver::read(const dictionary& dict)
{
    lagrangianType_ =
        lagrangianTypesNames_
        [
            dict.lookupOrDefault<word>("lagrangianType", "basicReactingCloud")
        ];
    cloudDictName_ = dict.lookupOrDefault<word>("dictName", "reactingCloud1");
}


void Foam::fv::lagrangianSolver::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName && corrector == 0)
    {
        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                parcels_().evolve();
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                multiphseParcels_().evolve();
                break;
            }
            case BasicThermoCloud:
            {
                basicThermoCloud_().evolve();
                break;
            }
            case CoalCloud:
            {
                coalCloud_().evolve();
                break;
            }
            case BasicSprayCloud:
            {
                basicSprayCloud_().evolve();
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                basicKinematicCollidingCloud_().evolve();
                break;
            }
            case BasicKinematicCloud:
            {
                basicKinematicCloud_().evolve();
                break;
            }
        }
    }
}


bool Foam::fv::lagrangianSolver::initialise()
{
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    const volScalarField& rho = obr_.lookupObject<volScalarField>("rho");
    const volVectorField& U = obr_.lookupObject<volVectorField>("U");
    g_ =
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                mesh_.time().constant(),
                mesh_.time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

    Info<< "Creating parces" << endl;
    switch (lagrangianType_)
    {
        case BasicReactingCloud:
        {
            parcels_ =
                new basicReactingCloud
                (
                    cloudDictName_,
                    rho,
                    U,
                    g_(),
                    slgThermo()
                );
            break;
        }
        case BasicReactingMultiphaseCloud:
        {
            multiphseParcels_ =
                new basicReactingMultiphaseCloud
                (
                    cloudDictName_,
                    rho,
                    U,
                    g_(),
                    slgThermo()
                );
            break;
        }
        case BasicThermoCloud:
        {
            basicThermoCloud_ =
                new basicThermoCloud
                (
                    cloudDictName_,
                    rho,
                    U,
                    g_(),
                    slgThermo()
                );
            break;
        }
        case CoalCloud:
        {
            coalCloud_ =
                new coalCloud(cloudDictName_, rho, U, g_(), slgThermo());
            break;
        }
        case BasicSprayCloud:
        {
            basicSprayCloud_ =
                new basicSprayCloud(cloudDictName_, rho, U, g_(), slgThermo());
            break;
        }
        case BasicKinematicCollidingCloud:
        {
            basicKinematicCollidingCloud_ =
                new basicKinematicCollidingCloud
                (
                    cloudDictName_,
                    rho,
                    U,
                    obr_.lookupObject<volScalarField>("mu"),
                    g_()
                );
            break;
        }
        case BasicKinematicCloud:
        {
            basicKinematicCloud_ =
                new basicKinematicCloud
                (
                    cloudDictName_,
                    rho,
                    U,
                    obr_.lookupObject<volScalarField>("mu"),
                    g_()
                );
            break;
        }
    }

    return true;
}


void Foam::fv::lagrangianSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {name()};

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::lagrangianSolver::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (fieldI == 0)
    {
        volScalarField& rho = obr_.lookupObjectRef<volScalarField>("rho");

        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                eqn += parcels_().Srho(rho);
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                eqn += multiphseParcels_().Srho(rho);
                break;
            }
            case BasicThermoCloud:
            {
                break;
            }
            case CoalCloud:
            {
                eqn += coalCloud_().Srho(rho);
                break;
            }
            case BasicSprayCloud:
            {
                eqn += basicSprayCloud_().Srho(rho);
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                break;
            }
            case BasicKinematicCloud:
            {
                break;
            }
        }
    }
}


void Foam::fv::lagrangianSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField& U = obr_.lookupObjectRef<volVectorField>("U");
    if (eqn.psi().name() == U.name())
    {
        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                eqn += parcels_().SU(U);
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                eqn += multiphseParcels_().SU(U);
                break;
            }
            case BasicThermoCloud:
            {
                eqn += basicThermoCloud_().SU(U);
                break;
            }
            case CoalCloud:
            {
                eqn += coalCloud_().SU(U);
                break;
            }
            case BasicSprayCloud:
            {
                eqn += basicSprayCloud_().SU(U);
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                break;
            }
            case BasicKinematicCloud:
            {
                break;
            }
        }
    }
}


void Foam::fv::lagrangianSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (fieldI == 0)
    {
        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                eqn += parcels_().Srho();
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                eqn += multiphseParcels_().Srho();
                break;
            }
            case BasicThermoCloud:
            {
                break;
            }
            case CoalCloud:
            {
                eqn += coalCloud_().Srho();
                break;
            }
            case BasicSprayCloud:
            {
                eqn += basicSprayCloud_().Srho();
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                break;
            }
            case BasicKinematicCloud:
            {
                break;
            }
        }
    }
    else if (fieldI == 2)
    {
        volScalarField& he = const_cast<volScalarField&>(eqn.psi());

        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                eqn += parcels_().Sh(he);
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                eqn += multiphseParcels_().Sh(he);
                break;
            }
            case BasicThermoCloud:
            {
                eqn += basicThermoCloud_().Sh(he);
                break;
            }
            case CoalCloud:
            {
                eqn += coalCloud_().Sh(he);
                break;
            }
            case BasicSprayCloud:
            {
                eqn += basicSprayCloud_().Sh(he);
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                break;
            }
            case BasicKinematicCloud:
            {
                break;
            }
        }
    }
    else if (fieldI > 2)
    {
        volScalarField& Yi = const_cast<volScalarField&>(eqn.psi());
        switch (lagrangianType_)
        {
            case BasicReactingCloud:
            {
                eqn += parcels_().SYi(fieldI - 3, Yi);
                break;
            }
            case BasicReactingMultiphaseCloud:
            {
                eqn += multiphseParcels_().SYi(fieldI - 3, Yi);
                break;
            }
            case BasicThermoCloud:
            {
                break;
            }
            case CoalCloud:
            {
                eqn += coalCloud_().SYi(fieldI - 3, Yi);
                break;
            }
            case BasicSprayCloud:
            {
                eqn += basicSprayCloud_().SYi(fieldI - 3, Yi);
                break;
            }
            case BasicKinematicCollidingCloud:
            {
                break;
            }
            case BasicKinematicCloud:
            {
                break;
            }
        }
    }
}


void Foam::fv::lagrangianSolver::correct
(
    const word& solveName,
    const word& regionName
)
{}


void Foam::fv::lagrangianSolver::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    DynamicList<word> fn;
    fn.append("rho");
    fn.append("U");
    fn.append((*thermoPtr_).he().name());

    if
    (
        lagrangianType_ != BasicKinematicCollidingCloud
     && lagrangianType_ != BasicKinematicCloud
    )
    {
        auto& Y =
            refCast<fluidMulticomponentThermo>(*thermoPtr_).composition().Y();

        forAll(Y, Yi)
        {
            if
            (
                Yi
             != refCast<fluidMulticomponentThermo>
                (
                    *thermoPtr_
                ).composition().defaultSpecie()
            )
            {
                fn.append(Y[Yi].name());
            }
        }
    }
    fields.transfer(fn);

    forAll(fields, i)
    {
        sourceDependencies.insert(fields[i], {"lagrangian"});
    }
}


scalar Foam::fv::lagrangianSolver::getMaxTimeStep()
{
    // Needs to be checked
    return GREAT;
}


// ************************************************************************* //
