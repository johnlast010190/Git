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
    (c) 2004-6 H. Jasak All rights reserved
    (c) 2024 Engys Ltd.

Class
    SegregatedSolver

Description
    Segregated solver for block matrices.  Segregation decouples block
    coefficients and calls a scalar solver for each component separately.
    Cross-component coupling is lagged in rhs.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    SegregatedSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SegregatedSolver_H
#define SegregatedSolver_H

#include "matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/blockLduSolvers.H"
#include "matrices/blockLduMatrix/BlockLduSolvers/BlockIterativeSolver/BlockIterativeSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SegregatedSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SegregatedSolver
:
    public BlockLduSolver<Type>
{
    // Private data

        //- Scalar x
        scalarField scalarX_;

        //- Scalar matrix
        blockScalarMatrix scalarMatrix_;

        //- Scalar b
        scalarField scalarB_;


public:

    //- Runtime type information
    TypeName("Segregated");


    // Constructors

        //- Construct from matrix components and solver data stream
        SegregatedSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        SegregatedSolver(const SegregatedSolver<Type>&) = delete;


        // Destructor

            virtual ~SegregatedSolver()
            {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<Type> solve
        (
            Field<Type>& x,
            const Field<Type>& b
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SegregatedSolver<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "matrices/blockLduMatrix/BlockLduSolvers/Segregated/SegregatedSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
