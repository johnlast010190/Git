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
    (c) 2018-2023 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

Class
    Foam::freestreamVelocityFvPatchVectorField

Description
    This boundary condition provides a free-stream condition for velocity.

    It is an inlet-outlet condition that uses the velocity orientation to
    continuously blend between fixed value for normal inlet and zero gradient
    for normal outlet flow.

Usage
    \table
        Property        | Description             | Required    | Default value
        freestreamValue | freestream velocity     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            freestreamVelocity;
        freestreamValue uniform (300 0 0);
    }
    \endverbatim

    Note:
        This condition is designed to operate with the freestreamPressure
        condition.

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    freestreamVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef freestreamVelocityFvPatchVectorField_H
#define freestreamVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class freestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class freestreamVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    autoPtr<Function1<vector>> freestreamValue_;

public:

    //- Runtime type information
    TypeName("freestreamVelocity");


    // Constructors

        //- Construct from patch and internal field
        freestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given freestreamVelocityFvPatchVectorField onto
        //  a new patch
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&
        );

        //- Copy constructor setting internal field reference
        freestreamVelocityFvPatchVectorField
        (
            const freestreamVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchVectorField>
            (
                new freestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

            const vectorField& freestreamValue() const
            {
                return refValue();
            }

            vectorField& freestreamValue()
            {
                return refValue();
            }

            //- Used within constrain HbyA, significantly
            //  reduces p-U decoupling
            virtual bool assignable() const override
            {
                return true;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            ) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&) override;

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;

            //- Returns frame velocity - overide and return zero
            virtual tmp<Field<vector>> getFrameVelocity() const override;

            //- Make relative doesn't do anything for this boundary condition
            virtual void makeRelative(scalarField&) const override
            {}

            //- Make relative doesn't do anything for this boundary condition
            virtual void makeRelative(vectorField&) const override
            {}


        //- Write
        virtual void write(Ostream&) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
