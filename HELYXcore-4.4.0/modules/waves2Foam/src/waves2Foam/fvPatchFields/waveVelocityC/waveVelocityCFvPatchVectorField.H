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
    (c) held by original author

Class
    Foam::waveVelocityCFvPatchVectorField

Description
    Boundary condition for the velocity field in a multiphase solver.
    To be used for wave boundary conditions and is catered towards ship capsize
    simulations where the domain is rotating significantly. Boundary condition
    uses a user-specified threshold to switch between an inlet and atmosphere
    boundary condition as the mesh moves.

    The boundary condition is derived from mixedFvPatchField and it is generic
    in the sense, that it loads the abstract class waveTheory. The boundary
    condition is based off the waveVelocity used in the original waves2Foam
    library.

SourceFiles
    waveVelocityCFvPatchVectorField.C

Author
    Derived from: Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.
    Created by: Kevin Silva, CSHL, The University of Michigan

\*---------------------------------------------------------------------------*/

#ifndef waveVelocityCFvPatchVectorField_H
#define waveVelocityCFvPatchVectorField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"
#include "convexPolyhedral/convexPolyhedral.H"
#include "waveTheories/waveTheory/waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class waveVelocityCFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveVelocityCFvPatchVectorField
:
    public mixedFvPatchField<vector>,
    public convexPolyhedral
{

private:

        // Private objects
        autoPtr<waveTheories::waveTheory> waveProps_;

        //- Flux field name
        word phiName_;

        //- Threshold for the dynamic boundary condition switch
        scalar thresholdValue;

        //- Threshold axis for the dynamic boundary condition switch (x,y,z)
        word thresholdAxis;

    // Private member functions

    //- Returns a scalarField of the signed distance to an arbitrary surface
        virtual void signedPointToSurfaceDistance
        (
            const pointField&,
            scalarField&
        );


    //- Returns a scalar of the signed distance to an arbitrary surface
        virtual scalar signedPointToSurfaceDistance
        (
            const point&
        ) const;

public:

    //- Runtime type information
    TypeName("waveVelocityC");


    // Constructors

        //- Construct from patch and internal field
        waveVelocityCFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveVelocityCFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given surfacewaveVelocityCFvPatchVectorField
        //  onto a new patch
        waveVelocityCFvPatchVectorField
        (
            const waveVelocityCFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new waveVelocityCFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveVelocityCFvPatchVectorField
        (
            const waveVelocityCFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector,
            volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new waveVelocityCFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "surfacewaveVelocityCFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
