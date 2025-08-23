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
    (c) 2022-2025 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "phaseChangeSource.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "binaryPhaseModels/phaseChangeModels/phaseChangeModel/phaseChangeModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseChangeSource, 0);
}
}

makeFvSolverOption(phaseChangeSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseChangeSource::phaseChangeSource
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

void Foam::fv::phaseChangeSource::read(const dictionary& dict)
{}

bool Foam::fv::phaseChangeSource::initialise()
{
    thermoPtr_ =
        dynamic_cast<multiphaseThermo*>(basicThermo::lookupPtr(obr_));

    if (!thermoPtr_)
    {
        Info<< "No phase change - not multiphase" << endl;
        return false;
    }

    const multiphaseThermo::phaseChangeModelTable& pcModels =
        thermoPtr_->phaseChangeModels();
    if (!pcModels.size())
    {
        Info<< "No phase change sources" << endl;
        return false;
    }
    else
    {
        Info<< "Phase change active" << endl;
    }

    return true;
}


void Foam::fv::phaseChangeSource::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("phaseChange");

    // Correct should be after mesh update
    optionalDependencies.insert(solveNames[0], {"fvMesh"});
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::phaseChangeSource::correct
(
    const word& solveName,
    const word& regionName
)
{
    forAllIters(thermoPtr_->phaseChangeModels(), iter)
    {
        iter()->correct();
    }
}


void Foam::fv::phaseChangeSource::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    fields = {"p"};

    const multiphaseThermo::phaseChangeModelTable& pcModels =
        thermoPtr_->phaseChangeModels();
    forAllConstIters(pcModels, iter)
    {
        // Add sources on both sides, but not passive phase
        const label passiveIndex = thermoPtr_->fractions().passiveIndex();
        if (thermoPtr_->phases().find(iter.key().first()) != passiveIndex)
        {
            fields.append(IOobject::groupName("alpha", iter.key().first()));
        }
        if (thermoPtr_->phases().find(iter.key().second()) != passiveIndex)
        {
            fields.append(IOobject::groupName("alpha", iter.key().second()));
        }
    }

    for (const word& f : fields)
    {
        sourceDependencies.insert(f, {"phaseChange"});
    }
}


void Foam::fv::phaseChangeSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // alpha equations
    const word thisPhaseName = eqn.psi().group();
    for (const word& otherPhaseName : thermoPtr_->phases())
    {
        if (otherPhaseName == thisPhaseName)
        {
            continue;
        }

        // Add source to alpha equation RHS
        const phaseChangeModel* pcModelPtr =
            thermoPtr_->phaseChangeModels(thisPhaseName, otherPhaseName);

        if (pcModelPtr)
        {
            Pair<tmp<volScalarField>> vDotAlphal = pcModelPtr->vDotAlphal();
            const volScalarField& vDotcAlphal = vDotAlphal[0]();
            const volScalarField& vDotvAlphal = vDotAlphal[1]();
            const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);
            eqn += fvm::Sp(vDotvmcAlphal, eqn.psi()) + vDotcAlphal;
        }

        // Add reverse of other-side source if opposite phase change exists
        const phaseChangeModel* otherPcModelPtr =
            thermoPtr_->phaseChangeModels(otherPhaseName, thisPhaseName);

        if (otherPcModelPtr)
        {
            Pair<tmp<volScalarField>> vDotAlphal =
                otherPcModelPtr->vDotAlphal();
            const volScalarField& vDotcAlphal = vDotAlphal[0]();
            const volScalarField& vDotvAlphal = vDotAlphal[1]();
            const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);
            eqn +=
                fvm::Sp(vDotvmcAlphal, eqn.psi())
                // The sum of alphas below will just be equal to 1 in the
                // two-phase case
              - vDotvmcAlphal*
                (
                    eqn.psi()+thermoPtr_->fractions().alphas(otherPhaseName)
                )
              - vDotcAlphal;
        }
    }
}


void Foam::fv::phaseChangeSource::addSup
(
    const volScalarField& psiByRho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    for (const word& phaseName : thermoPtr_->phases())
    {
        for (const word& otherPhaseName : thermoPtr_->phases())
        {
            if (otherPhaseName == phaseName)
            {
                continue;
            }

            const phaseChangeModel* pcModelPtr =
                thermoPtr_->phaseChangeModels
                (
                    phaseName, otherPhaseName
                );

            if (pcModelPtr)
            {
                // Add source to pressure equation RHS
                Pair<tmp<volScalarField>> vDotP = pcModelPtr->vDotP();
                const volScalarField& vDotcP = vDotP[0]();
                const volScalarField& vDotvP = vDotP[1]();

                eqn -=
                    fvm::Sp(vDotvP - vDotcP, thermoPtr_->p())
                  - (vDotvP - vDotcP)*pcModelPtr->pSat();
            }
        }
    }
}


// ************************************************************************* //
