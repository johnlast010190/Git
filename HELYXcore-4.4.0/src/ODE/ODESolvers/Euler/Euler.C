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
    (c) 2013-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "ODESolvers/Euler/Euler.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Euler, 0);
    addToRunTimeSelectionTable(ODESolver, Euler, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Euler::Euler(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    err_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Euler::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(err_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::Euler::solve
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    // Calculate error estimate from the change in state:
    forAll(err_, i)
    {
        err_[i] = dx*dydx0[i];
    }

    // Update the state
    forAll(y, i)
    {
        y[i] = y0[i] + err_[i];
    }

    return normalizeError(y0, y, err_);
}


void Foam::Euler::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    scalar& dxTry
) const
{
    adaptiveSolver::solve(odes_, x, y, li, dxTry);
}


// ************************************************************************* //
