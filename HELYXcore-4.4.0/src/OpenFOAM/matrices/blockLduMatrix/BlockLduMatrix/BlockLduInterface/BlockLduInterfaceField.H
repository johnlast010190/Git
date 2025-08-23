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
    (c) 2004-2006 Ivor Clifford
    (c) 2024 Engys Ltd.

Class
    BlockLduInterfaceField

Description
    Abstract base class for interface fields with block coefficients

Author
    Ivor Clifford, ivor.clifford@gmail.com

SourceFiles
    BlockLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterfaceField_H
#define BlockLduInterfaceField_H

#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
class CoeffField;

/*---------------------------------------------------------------------------*\
                           Class BlockLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterfaceField
:
    public lduInterfaceField
{
public:

    //- Runtime Type information
    TypeName("BlockLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        BlockLduInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}

        //- Disallow default bitwise copy construction
        BlockLduInterfaceField(const BlockLduInterfaceField&) = delete;


    // Destructor

        virtual ~BlockLduInterfaceField()
        {}



    // Member Functions

        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
