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
    (c) 2011 Engys Ltd.

Class
    unidirectionalVelocityOutletFvPatchVectorField

Description

SourceFiles
    unidirectionalVelocityOutletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef unidirectionalVelocityOutletFvPatchVectorField_H
#define unidirectionalVelocityOutletFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class unidirectionalVelocityOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class unidirectionalVelocityOutletFvPatchVectorField
:
    public inletOutletFvPatchVectorField
{
private:

    //- sink resistance factor for inflow next to the patch
    scalar sink_;

public:

    //- Runtime type information
    TypeName("unidirectionalVelocityOutlet");


    // Constructors

        //- Construct from patch and internal field
        unidirectionalVelocityOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unidirectionalVelocityOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given unidirectionalVelocityOutletFvPatchVectorField
        //  onto a new patch
        unidirectionalVelocityOutletFvPatchVectorField
        (
            const unidirectionalVelocityOutletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new unidirectionalVelocityOutletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unidirectionalVelocityOutletFvPatchVectorField
        (
            const unidirectionalVelocityOutletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new unidirectionalVelocityOutletFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //-Manipulate a matrix
        virtual void manipulateMatrix(fvVectorMatrix& matrix);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
