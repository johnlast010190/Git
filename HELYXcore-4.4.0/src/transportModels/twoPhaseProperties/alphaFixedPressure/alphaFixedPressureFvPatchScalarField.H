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

Class
    Foam::alphaFixedPressureFvPatchScalarField

Description
    A fixed-pressure alphaContactAngle boundary

SourceFiles
    alphaFixedPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaFixedPressureFvPatchScalarField_H
#define alphaFixedPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class alphaFixedPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaFixedPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Fixed pressure
        scalarField p_;


public:

    //- Runtime type information
    TypeName("alphaFixedPressure");


    // Constructors

        //- Construct from patch and internal field
        alphaFixedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaFixedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaFixedPressureFvPatchScalarField
        //  onto a new patch
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaFixedPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaFixedPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the alphaFixed pressure
            const scalarField& p() const
            {
                return p_;
            }

            //- Return reference to the alphaFixed pressure to allow adjustment
            scalarField& p()
            {
                return p_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
