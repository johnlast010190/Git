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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::offsetGAMGInterfaceField

Description
    GAMG agglomerated interface field.

SourceFiles
    offsetGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef offsetGAMGInterfaceField_H
#define offsetGAMGInterfaceField_H

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/GAMGInterfaceField/GAMGInterfaceField.H"
#include "offsetGAMGInterface.H"
#include "offsetLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class offsetGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class offsetGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public offsetLduInterfaceFieldBase
{

    // Private data members

        autoPtr< GAMGInterfaceField > slaveGAMGInterfaceField_;


public:

    //- Runtime type information
    TypeName("offset");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        offsetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        offsetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        offsetGAMGInterfaceField(const offsetGAMGInterfaceField&) = delete;


    //- Destructor
    virtual ~offsetGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return the underlying interface field
            virtual const lduInterfaceField& slaveInterfaceField() const
            {
                return slaveGAMGInterfaceField_();
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
                return slaveGAMGInterfaceField_->ready();
            }



        // Interface matrix update

            //- Initialise update based on interface functionality
            virtual void initInterfaceMatrixUpdate
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const offsetGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
