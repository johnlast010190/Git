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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::mappedFieldFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides a self-contained version of the \c mapped
    condition.  It does not use information on the patch; instead it holds
    the data locally.

Usage
    \table
        Property     | Description             | Required    | Default value
        field        | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedField;
        field           T;              // optional field name
        setAverage      no;             // apply an average value
        average         0;              // average to apply if setAverage
        value           uniform 0;      // place holder
    }
    \endverbatim

Note
    Since this condition can be applied on a per-field and per-patch basis,
    it is possible to duplicate the mapping information.  If possible, employ
    the \c mapped condition in preference to avoid this situation, and only
    employ this condition if it is not possible to change the underlying
    geometric (poly) patch type to \c mapped.

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchField

SourceFiles
    mappedFieldFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFieldFvPatchField_H
#define mappedFieldFvPatchField_H

#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "fields/fvPatchFields/derived/mappedField/mappedPatchFieldBase.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class mappedFieldFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFieldFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public mappedPatchBase,
    public mappedPatchFieldBase<Type>
{
protected:

    // Protected Member Functions

        virtual scalar getUfromFlowRate() const override;

public:

    //- Runtime type information
    TypeName("mappedField");


    // Constructors

        //- Construct from patch and internal field
        mappedFieldFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFieldFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct from patch, internal field and distance for normal type
        //  sampling
        mappedFieldFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,

            // mappedPatchBase
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance,

            // My settings
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct by mapping given
        //  mappedFieldFvPatchField
        //  onto a new patch
        mappedFieldFvPatchField
        (
            const mappedFieldFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedFieldFvPatchField
        (
            const mappedFieldFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        mappedFieldFvPatchField
        (
            const mappedFieldFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&) override;


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;


        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
