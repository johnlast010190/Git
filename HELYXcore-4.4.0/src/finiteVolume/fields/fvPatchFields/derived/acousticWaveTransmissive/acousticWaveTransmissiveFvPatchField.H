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
    Foam::acousticWaveTransmissiveFvPatchField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a wave transmissive outflow condition,
    based on solving DDt(W, field) = 0 at the boundary \c W is the wave velocity
    and \c field is the field to which this boundary condition is applied.

Usage
    \table
        Property       | Description             | Required    | Default value
        advectiveSpeed | Advective speed         | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            acousticWaveTransmissive;
        advectiveSpeed  343;
    }
    \endverbatim

See also
    Foam::advectiveFvPatchField

SourceFiles
    acousticWaveTransmissiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef acousticWaveTransmissiveFvPatchField_H
#define acousticWaveTransmissiveFvPatchField_H

#include "fields/fvPatchFields/derived/advective/advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class acousticWaveTransmissiveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class acousticWaveTransmissiveFvPatchField
:
    public advectiveFvPatchField<Type>
{

    // Private data

        //- Advection speed value
        scalar advectiveU_;


public:

    //- Runtime type information
    TypeName("acousticWaveTransmissive");


    // Constructors

        //- Construct from patch and internal field
        acousticWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        acousticWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given acousticWaveTransmissiveFvPatchField
        //  onto a new patch
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new acousticWaveTransmissiveFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new acousticWaveTransmissiveFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;

            //- Update the coefficients associated with the patch field
            //virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
