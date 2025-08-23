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
    BlockLduPrecon

Description
    Block LDU matrix preconditioner virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    newBlockLduPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduPrecon_H
#define BlockLduPrecon_H

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockLduPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduPrecon
{
protected:

    // Protected data

        //- Matrix reference
        const BlockLduMatrix<Type>& matrix_;


    // Protected member functions

        //- Find the smoother name (directly or from a sub-dictionary)
        static word getName(const dictionary&);


public:

    //- Runtime type information
    TypeName("BlockLduPrecon");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduPrecon,
            dictionary,
            (
                const BlockLduMatrix<Type>& matrix,
                const dictionary& dict,
                const word fieldName
            ),
            (
                matrix,
                dict,
                fieldName
            )
        );


    // Constructors

        //- Construct from matrix
        BlockLduPrecon(const BlockLduMatrix<Type>& matrix)
        :
            matrix_(matrix)
        {}

        //- Disallow default bitwise copy construction
        BlockLduPrecon(const BlockLduPrecon&) = delete;


    // Selectors

        //- Select given matrix and dictionary
        static autoPtr<BlockLduPrecon<Type>> New
        (
            const BlockLduMatrix<Type>&,
            const dictionary&,
            const word
        );


    // Destructor

        virtual ~BlockLduPrecon()
        {}


    // Member Functions

        const BlockLduMatrix<Type>& matrix() const
        {
            return matrix_;
        }

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const = 0;

        //- Execute preconditioning on a transposed matrix
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlockLduPrecon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "matrices/blockLduMatrix/BlockLduPrecons/BlockLduPrecon/BlockLduPrecon.C"
#endif

#define makeBlockPrecon(PreconType, typePreconType)                          \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typePreconType, 0);                      \
                                                                             \
addToRunTimeSelectionTable(PreconType, typePreconType, dictionary);

#define makeBlockPrecons(preconType)                                         \
                                                                             \
makeBlockPrecon(blockScalarPrecon, preconType##Scalar);                      \
makeBlockPrecon(blockVectorPrecon, preconType##Vector);                      \
makeBlockPrecon(blockTensorPrecon, preconType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
