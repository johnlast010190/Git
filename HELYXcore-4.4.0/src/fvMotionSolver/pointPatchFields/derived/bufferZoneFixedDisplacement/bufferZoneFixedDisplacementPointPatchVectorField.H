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
    (c) 2013-2024 Engys Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::bufferZoneFixedDisplacementPointPatchVectorField

Description
    Defines a rigid body rotation within a buffer zone a specified
    distance from the wall

SourceFiles
    bufferZoneFixedDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bufferZoneFixedDisplacementPointPatchVectorField_H
#define bufferZoneFixedDisplacementPointPatchVectorField_H

#include "pointPatchFields/derived/bufferZoneDisplacementBase/bufferZoneDisplacementBasePointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bufferZoneFixedDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bufferZoneFixedDisplacementPointPatchVectorField
:
    public bufferZoneDisplacementBasePointPatchVectorField
{
    // Private data

        vector radialVelocity_;
        vector linearVelocity_;
        vector CofG_;


public:

    //- Runtime type information
    TypeName("bufferZoneFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const bufferZoneFixedDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const override
        {
            return autoPtr<pointPatchField<vector>>
            (
                new bufferZoneFixedDisplacementPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const bufferZoneFixedDisplacementPointPatchVectorField&,
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
                new bufferZoneFixedDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

        // Evaluation functions

            //- Implement motion function
            virtual tmp<vectorField> getDisplacement
            (
                const vectorField& oldPoints
            ) override;

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
