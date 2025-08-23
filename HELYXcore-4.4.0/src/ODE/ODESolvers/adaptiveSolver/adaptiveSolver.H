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

Class
    Foam::adaptiveSolver

Group
    grpODESolvers

Description

SourceFiles
    adaptiveSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveSolver_H
#define adaptiveSolver_H

#include "ODESolvers/ODESolver/ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveSolver Declaration
\*---------------------------------------------------------------------------*/

class adaptiveSolver
{
    // Private data

        //- Step-size adjustment controls
        scalar safeScale_, alphaInc_, alphaDec_, minScale_, maxScale_;

        //- Cache for dydx at the initial time
        mutable scalarField dydx0_;

        //- Temprorary for the test-step solution
        mutable scalarField yTemp_;


public:

    // Constructors

        //- Construct from ODESystem
        adaptiveSolver(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~adaptiveSolver()
    {}


    // Member Functions

        //- Resize the ODE solver
        bool resize(const label n);

        //- Solve a single step dx and return the error
        virtual scalar solve
        (
            const scalar x0,
            const scalarField& y0,
            const label li,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const = 0;

        //- Solve the ODE system and the update the state
        void solve
        (
            const ODESystem& ode,
            scalar& x,
            scalarField& y,
            const label li,
            scalar& dxTry
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
