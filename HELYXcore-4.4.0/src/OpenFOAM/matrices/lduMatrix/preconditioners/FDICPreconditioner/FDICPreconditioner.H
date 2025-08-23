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
    Foam::FDICPreconditioner

Group
    grpLduMatrixPreconditioners

Description
    Faster version of the DICPreconditioner diagonal-based incomplete
    Cholesky preconditioner for symmetric matrices
    (symmetric equivalent of DILU) in which the the reciprocal of the
    preconditioned diagonal and the upper coefficients divided by the diagonal
    are calculated and stored.

SourceFiles
    FDICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef FDICPreconditioner_H
#define FDICPreconditioner_H

#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FDICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class FDICPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- The reciprocal preconditioned diagonal
        scalarField rD_;
        scalarField rDuUpper_;
        scalarField rDlUpper_;


public:

    //- Runtime type information
    TypeName("FDIC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        FDICPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );

        //- Disallow default bitwise copy construction
        FDICPreconditioner(const FDICPreconditioner&) = delete;


    //- Destructor
    virtual ~FDICPreconditioner()
    {}


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt=0
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FDICPreconditioner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
