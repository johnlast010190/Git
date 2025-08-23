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
    (c) 2013 OpenFOAM Foundation

Class
    Foam::interstitialInletVelocityFvPatchVectorField

Description
    Inlet velocity in which the actual interstitial velocity is calculated
    by dividing the specified inletVelocity field with the local phase-fraction.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type              interstitialInletVelocity;
        inletVelocity     uniform (0 0.2 0);// Non-interstitial inlet velocity
        alpha             alpha.particles;  // Name of the phase-fraction field
        value             uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    interstitialInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef interstitialInletVelocityFvPatchVectorField_H
#define interstitialInletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class interstitialInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class interstitialInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet velocity
        vectorField inletVelocity_;

        //- Name of the flux transporting the field
        word alphaName_;


public:

   //- Runtime type information
   TypeName("interstitialInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        interstitialInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interstitialInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interstitialInletVelocityFvPatchVectorField
        //  onto a new patch
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new interstitialInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&,
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
                new interstitialInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


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
