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
    (c) 2021 Engys Ltd.
    (c) 2017-2019, The Regents of the University of California
    (c) 2011-2018 OpenFOAM Foundation

Class
    turbulentATSMABLInletFvPatchVectorField

SourceFiles
    turbulentATSMABLInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentATSMABLInletFvPatchVectorField_H
#define turbulentATSMABLInletFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/derived/turbulentATSMInlet/turbulentATSMInletFvPatchVectorField.H"
#include "cfdTools/general/ABLProfile/ABLProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace incompressible
// {

/*---------------------------------------------------------------------------*\
       Class turbulentATSMABLInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentATSMABLInletFvPatchVectorField
:
    public turbulentATSMInletFvPatchVectorField
{

private:

    // Private data

        //- profile mode RichardsHoxey, AIJ, userDefined
        autoPtr<ABLProfile> profile_;

    // Private member functions

        //- initialise base class member data
        void initialiseBaseClassData
        (
            const fvPatch&, const dictionary&
        );

public:

    //- Runtime type information
    TypeName("turbulentATSMABLInlet");


    // Constructors

            //- Construct from patch and internal field
        turbulentATSMABLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentATSMABLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentATSMABLInletFvPatchVectorField
        //- onto a new patch
        turbulentATSMABLInletFvPatchVectorField
        (
            const turbulentATSMABLInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentATSMABLInletFvPatchVectorField
        (
            const turbulentATSMABLInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentATSMABLInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentATSMABLInletFvPatchVectorField
        (
            const turbulentATSMABLInletFvPatchVectorField&,
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
                new turbulentATSMABLInletFvPatchVectorField(*this, iF)
            );
        }

    // Member functions


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
