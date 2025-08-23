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
    (c) 2004-2011 H. Jasak
    (c) 2024 Engys Ltd.

Class
    BlockLduInterface

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterface_H
#define BlockLduInterface_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/FieldFields/FieldField/FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

/*---------------------------------------------------------------------------*\
                     Class BlockLduInterface Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterface
{
public:

    //- Runtime type information
    TypeName("BlockLduInterface");


    // Constructors

        //- Disallow default bitwise copy construction
        BlockLduInterface(const BlockLduInterface&) = delete;


    // Destructor

        virtual ~BlockLduInterface()
        {}


    // Member Functions

        // Coupled interface matrix update

            //- Initialise matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const BlockLduMatrix<Type>& matrix,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const
            {}

            //- Update result based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const BlockLduMatrix<Type>& matrix,
                const FieldField<CoeffField, Type>& coeffs,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
