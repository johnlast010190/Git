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
    (c) 2013-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "reactionsSource.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "basicThermo/basicThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(reactionsSource, 0);
}
}

makeFvSolverOption(reactionsSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::reactionsSource::reactionsSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    isReactingFoam_(dict.lookupOrDefault<Switch>("isReactingFoam", false))
{
    Info<< "Creating reaction model\n" << endl;

    const compressible::turbulenceModel& turbModel =
        obr.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    reaction_ =
        combustionModel::New
        (
            refCast<fluidMulticomponentThermo>
            (
                multiphaseThermo::lookupOrCreate(obr, phaseName_)
            ),
            turbModel
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::reactionsSource::read(const dictionary& dict)
{
    isReactingFoam_ = dict.lookupOrDefault<Switch>("isReactingFoam", false);
}


bool Foam::fv::reactionsSource::initialise()
{
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);
    return true;
}


void Foam::fv::reactionsSource::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {"reactions"};

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::reactionsSource::correct
(
    const word& solveName,
    const word& regionName
)
{
    reaction_->correct();
}


void Foam::fv::reactionsSource::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    DynamicList<word> fn;
    fn.append(reaction_->thermo().he().name());
    auto& Y = reaction_->thermo().composition().Y();

    forAll(Y, Yi)
    {
        if (Yi != reaction_->thermo().composition().defaultSpecie())
        {
            fn.append(Y[Yi].name());
        }
    }
    fields.transfer(fn);

    forAll(fields, i)
    {
        sourceDependencies.insert(fields[i], {"reactions"});
    }
}


void Foam::fv::reactionsSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    bool superficial =
        reaction_->thermo().properties().lookupOrDefault
        (
            "superficialReactionRate",
            false
        );
    scalar eps = 1;
    if (superficial)
    {
        eps =
            reaction_->thermo().properties().lookupOrDefault<scalar>
            (
                "porosity",
                1
            );
    }

    if (fieldI > 0)
    {
        volScalarField& Yi = const_cast<volScalarField&>(eqn.psi());
        if (superficial)
        {
            eqn += 1/eps*reaction_->R(Yi);
        }
        else
        {
            eqn += reaction_->R(Yi);
        }
    }
    else
    {
        if (superficial)
        {
            eqn += reaction_->Qdot()/eps;
        }
        else
        {
            eqn += reaction_->Qdot();
        }
    }
}


bool Foam::fv::reactionsSource::setRDeltaT()
{
    if (fv::localEulerDdt::enabled(mesh_))
    {
        const dictionary& solnDict = solnControlPtr_->dict();

        // The local time step will be registred under mesh object registry
        volScalarField& rDeltaT = lookupOrConstructRDeltaT();

        // Maximum change in cell temperature per iteration
        // (relative to previous value)
        const scalar alphaTemp = solnDict.lookupOrDefault("alphaTemp", 0.05);

        // Heat release rate time scale
        if (alphaTemp < 1)
        {
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>("rho");
            const volScalarField& T = reaction_->thermo().T();
            tmp<volScalarField::Internal> rDeltaTT
            (
                volScalarField::Internal::New
                (
                    "Qdot-rDeltaT",
                    mag(reaction_->Qdot())
                   /(alphaTemp*rho*reaction_->thermo().Cp()*T)
                )
            );

            Info<< "    Heat release rate = "
                << 1/(gMax(rDeltaTT().field()) + VSMALL) << ", "
                << 1/(gMin(rDeltaTT().field()) + VSMALL) << endl;

            rDeltaT.ref() = max(rDeltaT(), rDeltaTT);
        }

        // Combustion rate time scale

        // Maximum change in cell concentration per iteration
        // (relative to reference value)
        const scalar alphaY = solnDict.lookupOrDefault("alphaY", 1.0);
        if (alphaY < 1)
        {
            const dictionary Yref(solnDict.subDict("Yref"));

            tmp<volScalarField::Internal> rDeltaTY
            (
                volScalarField::Internal::New
                (
                    "rDeltaTY",
                    mesh_,
                    dimensionedScalar(rDeltaT.dimensions(), 0)
                )
            );

            bool foundY = false;
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>("rho");
            speciesMassFractions& composition =
                const_cast<speciesMassFractions&>
                (
                    dynamic_cast<const speciesMassFractions&>
                    (
                        reaction_->thermo().composition()
                    )
                );
            PtrList<volScalarField>& Y = composition.Y();

            const word speciesConcentrationName
            (
                obr_.names<speciesConcentrationSolver>().first()
            );
            const speciesConcentrationSolver* concentrationSolver =
                obr_.lookupObjectPtr<speciesConcentrationSolver>
                (
                    speciesConcentrationName
                );

            const label inertIndex = concentrationSolver->inertIndex();

            forAll(Y, i)
            {
                if (i != inertIndex && composition.active(i))
                {
                    volScalarField& Yi = Y[i];

                    if (Yref.found(Yi.name()))
                    {
                        foundY = true;
                        const scalar Yrefi = Yref.lookup<scalar>(Yi.name());

                        rDeltaTY.ref().field() =
                            max
                            (
                                mag
                                (
                                    reaction_->R(Yi)().source()
                                   /((Yrefi*alphaY)*(rho*mesh_.V()))
                                ),
                                rDeltaTY()
                            );
                    }
                }
            }

            if (foundY)
            {
                Info<< "    Species reactions = "
                    << 1/(gMax(rDeltaTY().field()) + VSMALL) << ", "
                    << 1/(gMin(rDeltaTY().field()) + VSMALL) << endl;

                rDeltaT.ref() = max(rDeltaT(), rDeltaTY);
            }
        }
        return true;
    }

    return false;
}


// ************************************************************************* //
