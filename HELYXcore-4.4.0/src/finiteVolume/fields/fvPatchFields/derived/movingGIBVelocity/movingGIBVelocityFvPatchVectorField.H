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
    (c) 2018 Engys Ltd.

Class
    Foam::movingGIBVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition is used on the GIB interface where
    the dynamic mesh calculates the boundary velocity specific to application
    (burning interface etc)

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            movingGIBVelocity;
        value           uniform (0 0 0);
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    movingGIBVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingGIBVelocityFvPatchVectorField_H
#define movingGIBVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingGIBVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingGIBVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

       //- Name of velocity field
       word UName_;


public:

    //- Runtime type information
    TypeName("movingGIBVelocity");


    // Constructors

        //- Construct from patch and internal field
        movingGIBVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingGIBVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingGIBVelocityFvPatchVectorField
        //  onto a new patch
        movingGIBVelocityFvPatchVectorField
        (
            const movingGIBVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingGIBVelocityFvPatchVectorField
        (
            const movingGIBVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingGIBVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingGIBVelocityFvPatchVectorField
        (
            const movingGIBVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new movingGIBVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
