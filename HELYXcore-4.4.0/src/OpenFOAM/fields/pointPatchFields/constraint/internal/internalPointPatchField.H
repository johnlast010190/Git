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
    (c) 2021 OpenFOAM Foundation

Class
    Foam::internalPointPatchField

Description
    Constraint pointPatchField to hold values for internal face exposed by
    sub-setting.

SourceFiles
    internalPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef internalPointPatchField_H
#define internalPointPatchField_H

#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "meshes/pointMesh/pointPatches/constraint/internal/internalPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class internalPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class internalPointPatchField
:
    public pointPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(internalPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        internalPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        internalPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        internalPointPatchField
        (
            const internalPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        internalPointPatchField
        (
            const internalPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new internalPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new internalPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return the constraint type this pointPatchField implements
        virtual const word& constraintType() const override
        {

            return internalPointPatch::typeName;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
