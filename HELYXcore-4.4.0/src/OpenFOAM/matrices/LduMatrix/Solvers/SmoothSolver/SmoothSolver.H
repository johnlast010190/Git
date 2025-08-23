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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::SmoothSolver

Description
    Iterative solver for symmetric and assymetric matrices which uses a
    run-time selected smoother e.g. GaussSeidel to converge the solution to
    the required tolerance.  To improve efficiency, the residual is evaluated
    after every nSweeps smoothing iterations.

SourceFiles
    SmoothSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SmoothSolver_H
#define SmoothSolver_H

#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SmoothSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class SmoothSolver
:
    public LduMatrix<Type, DType, LUType>::solver
{

protected:

    // Protected data

        //- Number of sweeps before the evaluation of residual
        label nSweeps_;

        //- Read the control parameters from the controlDict_
        virtual void readControls();


public:

    //- Runtime type information
    TypeName("SmoothSolver");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        SmoothSolver
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );


    // Member Functions

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/LduMatrix/Solvers/SmoothSolver/SmoothSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
