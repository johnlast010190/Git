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
    (c) 2024 Engys Ltd.

Class
    Foam::PBiCCCG

Description
    Preconditioned bi-conjugate gradient solver for asymmetric lduMatrices
    using a run-time selectable preconditioner.

SourceFiles
    PBiCCCG.C

\*---------------------------------------------------------------------------*/

#ifndef PBiCCCG_H
#define PBiCCCG_H

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PBiCCCG Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class PBiCCCG
:
    public LduMatrix<Type, DType, LUType>::solver
{
public:

    //- Runtime type information
    TypeName("PBiCCCG");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        PBiCCCG
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );

        //- Disallow default bitwise copy construction
        PBiCCCG(const PBiCCCG&) = delete;


    // Destructor

        virtual ~PBiCCCG()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PBiCCCG&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/LduMatrix/Solvers/PBiCCCG/PBiCCCG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
