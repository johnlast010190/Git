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
    Foam::sedimentationPointPatchVectorField

Description
    Foam::sedimentationPointPatchVectorField

SourceFiles
    sedimentationPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sedimentationPointPatchVectorField_H
#define sedimentationPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class sedimentationPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sedimentationPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
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

        //- Is displacement?
        Switch isDisplacement_;


public:

    //- Runtime type information
    TypeName("sedimentation");


    // Constructors

        //- Construct from patch and internal field
        sedimentationPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sedimentationPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        sedimentationPointPatchVectorField
        (
            const sedimentationPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const override
        {
            return autoPtr<pointPatchField<vector>>
            (
                new sedimentationPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sedimentationPointPatchVectorField
        (
            const sedimentationPointPatchVectorField&,
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
                new sedimentationPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
