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

#include "pcweAcousticSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "algorithms/subCycle/subCycle.H"
#include "transport/speedSoundGas/speedSoundGas.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(pcweAcousticSolver, 0);
}
}

makeFvSolverOption(pcweAcousticSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::pcweAcousticSolver::pcweAcousticSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    thermoPtr_(nullptr),
    phiPtr_(nullptr),
    UPtr_(nullptr),
    pPtr_(nullptr),
    startTime_(0.0),
    dtAco_(mesh().time().deltaTValue()),
    psiAco_
    (
        volScalarField
        (
            IOobject
            (
                "psiAco",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    pAco_
    (
        volScalarField
        (
            IOobject
            (
                "pAcoustic",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("pAcoustic", dimPressure, 0.0)
        )
    ),
    dpdt_
    (
        volScalarField
        (
            IOobject
            (
                "dpdt",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("dpdt", dimArea/pow3(dimTime), 0.0)
        )
    ),
    convectiveTerms_(false),
    flowSolverSource_(false),
    nAcoSubCycles_(0)
{
    convectiveTerms_ = dict.lookupOrDefault<bool>("convectiveTerms", false);
    flowSolverSource_ = dict.lookupOrDefault<bool>("flowSolverSource", false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::pcweAcousticSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("psiAco");

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
    correctorMembers.insert
    (
        "nonOrthogonalCorrector:" + solveNames[0], solveNames
    );
}


void Foam::fv::pcweAcousticSolver::read(const dictionary& dict)
{
    startTime_ = dict.lookupOrDefault("startTime", 0.0);
    dtAco_ = dict.lookupOrDefault("dtAco", mesh().time().deltaTValue());
}


bool Foam::fv::pcweAcousticSolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    thermoPtr_ = &basicThermo::lookupOrCreate(obr_);
    c0_ = thermoPtr_->materials()(c0Model::typeName)();

    if (flowSolverSource_)
    {
        pPtr_ = &mesh_.lookupObject<volScalarField>("p");
    }

    if (convectiveTerms_)
    {
        phiPtr_ = &mesh_.lookupObjectRef<surfaceScalarField>("phi");
        UPtr_ = &mesh_.lookupObject<volVectorField>("U");
    }

    return true;
}


bool Foam::fv::pcweAcousticSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == "nonOrthogonalCorrector:" + psiAco_.name())
    {
        // Non-orthogonal correctors
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else if (correctorName == solverObject::outerCorrectorName)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::pcweAcousticSolver::solveEquation
(
    volScalarField& field
)
{
    fv::options& fvOptions = this->fvOptions();

    // Assemble Perturbed Convective Wave Equation
    tmp<fvScalarMatrix> tpcweEqn
    (
        new fvScalarMatrix
        (
            fvm::d2dt2(field)
           -sqr(c0_())*fvm::laplacian(field)
        ==
            fvOptions.d2dt2(field)
        )
    );

    if (convectiveTerms_)
    {
        const volVectorField& U(*UPtr_);

        tpcweEqn.ref() += (2*U & fvc::grad(fvc::ddt(field)))
            +((U & U)*fvm::laplacian(field));
    }

    tpcweEqn.ref() += dpdt_;

    tpcweEqn->relax();
    fvOptions.constrain(tpcweEqn.ref());

    tpcweEqn->solve();
    fvOptions.correct(field);

    if (convectiveTerms_)
    {
        pAco_ = thermoPtr_->rho()*fvc::ddt(field)
            + fvc::div(*phiPtr_, field);
    }
    else
    {
        pAco_ = thermoPtr_->rho()*fvc::ddt(field);
    }
}

void Foam::fv::pcweAcousticSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (time().value() > startTime_)
    {
        if (flowSolverSource_ && time().value() > startTime_)
        {
            const volScalarField& p(*pPtr_);
            dpdt_ = fvc::ddt(p)/thermoPtr_->rho();

            if (convectiveTerms_)
            {
                dpdt_ += fvc::div(*phiPtr_, p) / sqr(thermoPtr_->rho());
            }
        }

        nAcoSubCycles_ = nSubCycles();

        if (nAcoSubCycles_ > 1)
        {
            Time& time = const_cast<Time&>(mesh().time());

            {
                subCycleTime subCycle(time, nAcoSubCycles_);

                // Must reset time index to sub-cycling index due to
                // reset below
                psiAco_.timeIndex() = mesh_.time().timeIndex();

                const volScalarField source = (dpdt_/nAcoSubCycles_).ref();
                dpdt_ = source;

                while (!(++subCycle).end())
                {
                    solveEquation(psiAco_);
                    dpdt_ += source;
                }
            }

            // Protect old-times from modification outside the sub-cycles
            psiAco_.timeIndex() = mesh_.time().timeIndex();
        }
        else
        {
            solveEquation(psiAco_);
        }
    }
}


Foam::label Foam::fv::pcweAcousticSolver::nSubCycles() const
{
    const scalar deltaT = mesh_.time().deltaTValue();
    const label nSubCycles = round(deltaT/dtAco_);

    return nSubCycles;
}

// ************************************************************************* //
