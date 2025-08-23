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
    (c) 2022-2024 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "waves.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "relaxationZone/relaxationZone.H"
#include "cfdTools/general/include/fvCFD.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(waves, 0);
}
}

makeFvSolverOption(waves);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::waves::waves
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    thermoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::waves::read(const dictionary& dict)
{}

bool Foam::fv::waves::initialise()
{
    Info<< "Reading material properties\n" << endl;
    thermoPtr_ =
        &refCast<multiphaseThermo>(basicThermo::lookupOrCreate(obr_));

    if (thermoPtr_->alphas().size() != 2)
    {
        FatalErrorInFunction
            << "The " << this->typeName << " solver only supports two phase "
            << "flow." << nl << exit(FatalError);
    }
    alpha1_ =
        thermoPtr_->fractions().passiveIndex() == 0
      ? &thermoPtr_->alphas()[1]
      : &thermoPtr_->alphas()[0];

    relaxing_.set
    (
        new relaxationZone
        (
            mesh_, obr().lookupObjectRef<volVectorField>("U"), *alpha1_
        )
    );

    return true;
}


void Foam::fv::waves::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("waveRelaxation");

    // Relaxation correct should be after alpha solve and before fluid solves
    requiredDependencies.insert(solveNames[0], {alpha1_->name()});
    optionalDependencies.insert("U", solveNames);
    optionalDependencies.insert(thermoPtr_->p().name(), solveNames);

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::waves::correct
(
    const word& solveName,
    const word& regionName
)
{
    relaxing_().correct();
    thermoPtr_->correct();
}

// ************************************************************************* //
