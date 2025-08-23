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

#include "surfaceFilmSolver.H"
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
    defineTypeNameAndDebug(surfaceFilmSolver, 0);
}
}

makeFvSolverOption(surfaceFilmSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::surfaceFilmSolver::surfaceFilmSolver
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

void Foam::fv::surfaceFilmSolver::read(const dictionary& dict)
{}


bool Foam::fv::surfaceFilmSolver::initialise()
{
    if (!obr_.found(IOobject::groupName("U", phaseName_)))
    {
        Info<< "Reading field U\n" << endl;
        obr_.store
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", phaseName_),
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
            )
        );
    }
    if (!obr_.found(IOobject::groupName("rho", phaseName_)))
    {
        obr_.store
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho", phaseName_),
                    mesh_.time().timeName(),
                    obr_
                ),
                mesh_,
                dimensionedScalar(dimDensity, 0)
            )
        );
    }
    if (!obr_.found(IOobject::groupName("phi", phaseName_)))
    {
        obr_.store
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phi", phaseName_),
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT
                ),
                mesh_,
                dimensionedScalar(dimDensity, 0)
            )
        );
    }
    Info<< "\nConstructing surface film model" << endl;
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
    surfaceFilm_ =
        regionModels::surfaceFilmModel::New(mesh_, g_());
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    return true;
}


void Foam::fv::surfaceFilmSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {"surfaceFilm"};

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::surfaceFilmSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    surfaceFilm_->evolve();
}


void Foam::fv::surfaceFilmSolver::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    DynamicList<word> fn;
    fn.append("rho");
    fn.append((*thermoPtr_).he().name());
    auto& Y = refCast<fluidMulticomponentThermo>(*thermoPtr_).composition().Y();

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
    fields.transfer(fn);

    forAll(fields, i)
    {
        sourceDependencies.insert(fields[i], {"surfaceFilm"});
    }
}


void Foam::fv::surfaceFilmSolver::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (fieldI == 0)
    {
        eqn += surfaceFilm_().Srho();
    }
}


void Foam::fv::surfaceFilmSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (fieldI == 0)
    {
        eqn += surfaceFilm_().Srho();
    }
    else if (fieldI == 1)
    {
        eqn += surfaceFilm_().Sh();
    }
    else
    {
        eqn += surfaceFilm_().Srho(fieldI - 2);
    }
}


scalar Foam::fv::surfaceFilmSolver::getMaxTimeStep()
{
    // Needs to be implemented
    return GREAT;
}


// ************************************************************************* //
