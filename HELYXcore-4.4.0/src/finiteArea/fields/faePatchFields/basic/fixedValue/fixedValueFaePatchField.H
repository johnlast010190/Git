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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Engys Ltd.

Class
    fixedValueFaePatchField

Description

SourceFiles
    fixedValueFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueFaePatchField_H
#define fixedValueFaePatchField_H

#include "fields/faePatchFields/faePatchField/faePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedValueFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from patch and internal field
        fixedValueFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedValueFaePatchField<Type>
        //  onto a new patch
        fixedValueFaePatchField
        (
            const fixedValueFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueFaePatchField
        (
            const fixedValueFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>
            (
                new fixedValueFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedValueFaePatchField
        (
            const fixedValueFaePatchField<Type>&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, faEdgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type>>
            (
                new fixedValueFaePatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const faePatchField<Type>&) {}
        virtual void operator+=(const faePatchField<Type>&) {}
        virtual void operator-=(const faePatchField<Type>&) {}
        virtual void operator*=(const faePatchField<scalar>&) {}
        virtual void operator/=(const faePatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faePatchFields/basic/fixedValue/fixedValueFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
