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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2020-2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "helyxCoupledControl.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helyxCoupledControl, 0);

    addToRunTimeSelectionTable
    (
        solutionControl,
        helyxCoupledControl,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::helyxCoupledControl::readControls()
{
    solutionControl::readControls(false);

    const dictionary& solverDict = dict();

    solveFlow_ = solverDict.lookupOrDefault<Switch>("solveFlow", true);
    nCorr_ = solverDict.lookupOrDefault<label>("nOuterCorrectors", 1);
    turbOnFinalIterOnly_ =
        solverDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", false);
    stabilityMode_ = solverDict.lookupOrDefault<label>("stabilityMode", 1);
    choiCorrection_ =
        solverDict.lookupOrDefault<Switch>("choiCorrection", false);
}


bool Foam::helyxCoupledControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.blockSolverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi = applyToField(variableName);
        if (fieldi != -1)
        {

            scalar residual = 0;
            const scalar firstResidual =
                maxResidual(iter(), residual);

            checked = true;

            if (storeIni)
            {
                residualControl_[word::null][fieldi].initialResidual =
                    firstResidual;
            }

            const bool absCheck = residual < residualControl()[fieldi].absTol;
            bool relCheck = false;

            scalar relative = 0.0;
            if (!storeIni)
            {
                const scalar iniRes =
                    residualControl()[fieldi].initialResidual
                  + ROOTVSMALL;

                relative = residual/iniRes;
                relCheck = relative < residualControl()[fieldi].relTol;
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << variableName
                    << " outer iter " << corr_
                    << ": ini res = "
                    << residualControl()[fieldi].initialResidual
                    << ", abs tol = " << residual
                    << " (" << residualControl()[fieldi].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl()[fieldi].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::helyxCoupledControl::setFirstIterFlag
(
    const bool check,
    const bool force
)
{
    DebugInformation
        << "corr:" << corr_
        << endl;

    solutionControl::setFirstIterFlag(check && corr_ <= 1, force);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helyxCoupledControl::helyxCoupledControl
(
    fvMesh& mesh,
    const word& dictName
)
:
    helyxCoupledControl(mesh, mesh.thisDb(), dictName)
{
}


Foam::helyxCoupledControl::helyxCoupledControl
(
    fvMesh& mesh,
    const objectRegistry& obr,
    const word& dictName
)
:
    solutionControl(mesh, obr, dictName),
    solveFlow_(true),
    nCorr_(0),
    corr_(0),
    turbOnFinalIterOnly_(true),
    stabilityMode_(1),
    choiCorrection_(false),
    converged_(false)
{
    readControls();

    if (nCorr_ > 1)
    {
        Info<< nl;
        if (residualControl_.empty())
        {
            Info<< algorithmName_ << ": no residual control data found. "
                << "Calculations will employ " << nCorr_
                << " corrector loops" << nl << endl;
        }
        else
        {
            Info<< algorithmName_ << ": max iterations = " << nCorr_
                << endl;
            forAll(residualControl(), i)
            {
                Info<< "    field " << residualControl()[i].name << token::TAB
                    << ": relTol " << residualControl()[i].relTol
                    << ", tolerance " << residualControl()[i].absTol
                    << nl;
            }
            Info<< endl;
        }
    }
    else
    {
        Info<< nl << algorithmName_ << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helyxCoupledControl::~helyxCoupledControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::helyxCoupledControl::loop()
{
    readControls();

    corr_++;

    if (debug)
    {
        Info<< mesh_.time().timeName() << " loop: corr = " << corr_ << endl;
    }

    setFirstIterFlag();

    if (corr_ == nCorr_ + 1)
    {
        if ((!residualControl_.empty()) && (nCorr_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                 << nCorr_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                 << " iterations" << endl;

            mesh_.data::remove("finalIteration");

            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;

            mesh_.data::add("finalIteration", true);

            converged_ = true;
        }
    }
    else
    {
        if (corr_ <= nCorr_)
        {
            if (nCorr_!=1)
            {
                Info<< "Outer loop "
                     << mesh_.time().timeName()
                     << " | " << corr_ << endl;
            }

            if ((corr_ == nCorr_) && (nCorr_>1))
            {
                mesh_.data::add("finalIteration", true);
            }


            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
