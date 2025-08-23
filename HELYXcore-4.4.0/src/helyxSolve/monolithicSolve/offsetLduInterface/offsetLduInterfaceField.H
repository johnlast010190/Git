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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2015 Engys Ltd.

Class
    Foam::offsetLduInterfaceField

Description


SourceFiles
    offsetLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef offsetLduInterfaceField_H
#define offsetLduInterfaceField_H

#include "offsetLduInterface.H"
#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.H"
#include "offsetLduInterfaceFieldBase.H"
#include "fields/Fields/primitiveFieldsFwd.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class offsetLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class offsetLduInterfaceField
:
    public lduInterfaceField,
    virtual public offsetLduInterfaceFieldBase
{
    // Private data

        //- Underlying interface field that we are controlling
        const lduInterfaceField& slaveInterfaceField_;

        //- Offset into internal cell array of our internal field
        const label internalOffset_;

        //- Offset of the internal field coupled to (usually the same)
        const label foreignOffset_;

        //- Size of this region of internal cell array
        const label internalSize_;

        //- Size of coupled-to region of internal cell array
        const label foreignSize_;


public:

    //- Runtime type information
    TypeName("offset");


    // Constructors

        //- Construct from interface field and patch
        offsetLduInterfaceField
        (
            const lduInterfaceField& field,
            const lduInterface& slaveInterface,
            const lduInterface& offsetInterface,
            const label internalOffset,
            const label internalSize,
            const label foreignOffset,
            const label foreignSize
        );

        //- Disallow default bitwise copy construction
        offsetLduInterfaceField(const offsetLduInterfaceField&) = delete;


    //- Destructor
    virtual ~offsetLduInterfaceField();


    // Member Functions

        // Access

            //- Return the underlying interface field
            virtual const lduInterfaceField& slaveInterfaceField() const
            {
                return slaveInterfaceField_;
            }

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                //return slaveInterfaceField_.interfaceFieldType();
                return type();
            }


        // Coupled interface matrix update

            //- Is all data available
            virtual bool ready() const
            {
                return slaveInterfaceField_.ready();
            }


        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction dirn,
                const Pstream::commsTypes commsType
            ) const;


            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction dirn,
                const Pstream::commsTypes commsType
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const offsetLduInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
