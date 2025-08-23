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
    (c) 2024 Engys Ltd.

Description
    Sedimentation boundary condition for the drift scalar transport equation.

SourceFiles
    sedimentationFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef sedimentationFvPatchField_H
#define sedimentationFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sedimentationFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class sedimentationFvPatchField
:
    public fixedValueFvPatchField<scalar>
{
    // Private data

        //- Should the snow velocity/displacement be limited?
        Switch constrainDisplacement_;

        //- Name of the tri surface constrain
        word sedimentationSurfaceConstrainName_;

        //- Triangular mesh for the snow errosion limits
        autoPtr<triSurfaceMesh> sedimentationSurfaceConstrain_;

        //- Gravity direction
        vector normal_;

        //- Time
        scalar timeIndex_;

        //- Min value
        scalar minValue_;

        //- Max value
        scalar maxValue_;


public:

    //- Runtime type information
    TypeName("sedimentation");


    // Constructors

        //- Construct from patch and internal field
        sedimentationFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sedimentationFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        sedimentationFvPatchField
        (
            const sedimentationFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sedimentationFvPatchField
        (
            const sedimentationFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField > clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new sedimentationFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sedimentationFvPatchField
        (
            const sedimentationFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new sedimentationFvPatchField(*this, iF)
            );
        }


    //- Destructor
    virtual ~sedimentationFvPatchField();


    // Member functions

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

