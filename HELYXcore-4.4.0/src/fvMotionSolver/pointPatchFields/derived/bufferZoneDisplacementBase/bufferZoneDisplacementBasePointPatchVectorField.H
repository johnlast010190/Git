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
    Foam::bufferZoneDisplacementBasePointPatchVectorField

Description
    Base class for geometry morphing within a buffer zone defined though a
    specified distance from the wall.

SourceFiles
    bufferZoneDisplacementBasePointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bufferZoneDisplacementBasePointPatchVectorField_H
#define bufferZoneDisplacementBasePointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bufferZoneDisplacementBasePointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bufferZoneDisplacementBasePointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        pointField p0_;

        scalar bufferSize_;

        labelList bufferCells_;

        labelList bufferPoints_;


public:

    //- Runtime type information
    TypeName("bufferZoneDisplacementBase");


    // Constructors

        //- Construct from patch and internal field
        bufferZoneDisplacementBasePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bufferZoneDisplacementBasePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        bufferZoneDisplacementBasePointPatchVectorField
        (
            const bufferZoneDisplacementBasePointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const override
        {
            return autoPtr<pointPatchField<vector>>
            (
                new bufferZoneDisplacementBasePointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bufferZoneDisplacementBasePointPatchVectorField
        (
            const bufferZoneDisplacementBasePointPatchVectorField&,
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
                new bufferZoneDisplacementBasePointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const pointPatchFieldMapper&) override;

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            ) override;

        // Access functions

            //- Const access to p0_
            const pointField& p0()
            {
                return p0_;
            }

            //- Const access to bufferPoints_
            const labelList& bufferPoints()
            {
                return bufferPoints_;
            }

            //- Const access to bufferCells_
            const labelList& bufferCells()
            {
                return bufferCells_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;

            //- Set internal cells to constrain
            virtual void manipulateMatrix
            (
                const vectorField&,
                labelList&,
                vectorField&
            ) override;

            //- Set internal points to constrain
            virtual void setField
            (
                const vectorField&,
                vectorField&,
                bool setAll = false
            ) override;

            //- Return cells inside buffer region
            void calculateBufferCellsAndPoints();

            //- Get displacement given point field
            //  Motion descriptor for patch and internal points/cells
            virtual tmp<vectorField> getDisplacement
            (
                const vectorField& oldPoints
            );

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
