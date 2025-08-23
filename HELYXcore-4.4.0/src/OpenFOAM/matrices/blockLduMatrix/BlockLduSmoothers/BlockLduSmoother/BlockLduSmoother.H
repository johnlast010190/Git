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
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

SourceFiles
    newBlockLduSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduSmoother_H
#define BlockLduSmoother_H

#include "matrices/blockLduMatrix/BlockLduMatrix/blockLduMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BlockLduSmoother Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduSmoother
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
    TypeName("BlockLduSmoother");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockLduSmoother,
            dictionary,
            (
                const BlockLduMatrix<Type>& matrix,
                const dictionary& dict
            ),
            (
                matrix,
                dict
            )
        );


    // Constructors

        //- Construct from matrix
        BlockLduSmoother(const BlockLduMatrix<Type>& matrix)
        :
            matrix_(matrix)
        {}

        //- Disallow default bitwise copy construction
        BlockLduSmoother(const BlockLduSmoother&) = delete;


    // Selectors

        //- Select given matrix and dictionary
        static autoPtr<BlockLduSmoother<Type>> New
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockLduSmoother()
        {}


    // Member Functions

        //- Execute smoothing
        virtual void smoothCoarsest
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const
        {
            FatalErrorInFunction
                << "No smoothCoarsest, only valid from coarsest solvers"
                   ".\n"
                   "Something is wrong."
                << exit(FatalError);
        }

        //- Execute smoothing
        virtual void smooth
        (
            Field<Type>& x,
            const Field<Type>& b,
            const label nSweeps
        ) const = 0;

        virtual void initMatrix(){}

        //access matrix for debug
        const BlockLduMatrix<Type>& matrix() const
        {
            return matrix_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlockLduSmoother&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "matrices/blockLduMatrix/BlockLduSmoothers/BlockLduSmoother/BlockLduSmoother.C"
#endif

#define makeBlockSmoother(SmootherType, typeSmootherType)                    \
                                                                             \
defineNamedTemplateTypeNameAndDebug(typeSmootherType, 0);                    \
                                                                             \
addToRunTimeSelectionTable(SmootherType, typeSmootherType, dictionary);

#define makeBlockSmoothers(smootherType)                                       \
                                                                             \
makeBlockSmoother(blockScalarSmoother, smootherType##Scalar);                  \
makeBlockSmoother(blockVectorSmoother, smootherType##Vector);                  \
//makeBlockSmoother(blockTensorSmoother, smootherType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
