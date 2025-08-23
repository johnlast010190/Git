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

Class
    Foam::modalForcedMotionPointPatchVectorField

Description
    A forced motion implementation based on tabluar displacement input
    tailored to small rotor blade deformations.

SourceFiles
    modalForcedMotionPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef modalForcedMotionPointPatchVectorField_H
#define modalForcedMotionPointPatchVectorField_H

#include "pointPatchFields/derived/bufferZoneDisplacementBase/bufferZoneDisplacementBasePointPatchVectorField.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class modalForcedMotionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class modalForcedMotionPointPatchVectorField
:
    public bufferZoneDisplacementBasePointPatchVectorField
{
    // Private data

        //- Blade hub origin (centre of rotation)
        const vector origin_;

        //- Blade orientation
        const vector axis_;

        //- Radial displacement table (r, (dx,dy,dz))
        autoPtr<Function1<vector>> radialDisplacement_;

        //- Amplification of first mode (t, amp)
        autoPtr<Function1<scalar>> amplification_;


public:

    //- Runtime type information
    TypeName("modalForcedMotion");


    // Constructors

        //- Construct from patch and internal field
        modalForcedMotionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        modalForcedMotionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        modalForcedMotionPointPatchVectorField
        (
            const modalForcedMotionPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const override
        {
            return autoPtr<pointPatchField<vector>>
            (
                new modalForcedMotionPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        modalForcedMotionPointPatchVectorField
        (
            const modalForcedMotionPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const override
        {
            return autoPtr<pointPatchField<vector>>
            (
                new modalForcedMotionPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Compute radial coordinate for patch points
            tmp<scalarField> getRadialCoordinates
            (
                const vectorField& oldPoints
            );

            //- Compute max displacement field
            void calcMaxDisplacement
            (
                const vectorField& oldPoints,
                vectorField& mxDisplacement
            );

            //- Implement motion function
            virtual tmp<vectorField> getDisplacement
            (
                const vectorField& oldPoints
            ) override;

            // Return const access to amplification
            const Function1<scalar>& amplification() const
            {
                return amplification_;
            }

            // Return const access to amplification
            const Function1<vector>& radialDisplacement() const
            {
                return radialDisplacement_;
            }

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
