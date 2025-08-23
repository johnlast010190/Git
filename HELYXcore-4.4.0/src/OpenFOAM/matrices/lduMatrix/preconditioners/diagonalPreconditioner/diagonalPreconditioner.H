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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::diagonalPreconditioner

Group
    grpLduMatrixPreconditioners

Description
    Diagonal preconditioner for both symmetric and asymmetric matrices.

    The reciprocal of the diagonal is calculated and stored for reuse
    because on most systems '*' is faster than '/'.

SourceFiles
    diagonalPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef diagonalPreconditioner_H
#define diagonalPreconditioner_H

#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diagonalPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class diagonalPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- The reciprocal diagonal
        scalarField rD;


public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        diagonalPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );

        //- Disallow default bitwise copy construction
        diagonalPreconditioner(const diagonalPreconditioner&) = delete;


    //- Destructor
    virtual ~diagonalPreconditioner()
    {}


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt=0
        ) const;

        //- Return wT the transpose-matrix preconditioned form of residual rT.
        virtual void preconditionT
        (
            scalarField& wT,
            const scalarField& rT,
            const direction cmpt=0
        ) const
        {
            return precondition(wT, rT, cmpt);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const diagonalPreconditioner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
