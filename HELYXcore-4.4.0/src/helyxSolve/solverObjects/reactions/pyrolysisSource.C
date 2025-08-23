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

#include "pyrolysisSource.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "basicThermo/basicThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(pyrolysisSource, 0);
}
}

makeFvSolverOption(pyrolysisSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::pyrolysisSource::pyrolysisSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{
    Info<< "Creating pyrolysis model" << endl;
    pyrolysisPtr_ =
        new regionModels::pyrolysisModelCollection(mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::pyrolysisSource::read(const dictionary& dict)
{
    maxDi_ = pyrolysisPtr_->maxDiff();
    diNum_ = pyrolysisPtr_->solidRegionDiffNo();
}


bool Foam::fv::pyrolysisSource::initialise()
{
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    return true;
}


void Foam::fv::pyrolysisSource::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {"pyrolysis"};

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::pyrolysisSource::correct
(
    const word& solveName,
    const word& regionName
)
{
    maxDi_ = pyrolysisPtr_().maxDiff();
    diNum_ = pyrolysisPtr_().solidRegionDiffNo();
    pyrolysisPtr_().evolve();
}


void Foam::fv::pyrolysisSource::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{}


scalar Foam::fv::pyrolysisSource::getMaxTimeStep()
{
    return mesh_.time().deltaTValue()*min(maxDi_/(diNum_ + SMALL), 1.2);
}


// ************************************************************************* //
