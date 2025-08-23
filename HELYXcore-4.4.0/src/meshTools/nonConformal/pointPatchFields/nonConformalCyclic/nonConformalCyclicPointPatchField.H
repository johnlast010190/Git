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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

Class
    Foam::nonConformalCyclicPointPatchField

Description
    Constraint pointPatchField for nonConformalCyclic couplings.

SourceFiles
    nonConformalCyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicPointPatchField_H
#define nonConformalCyclicPointPatchField_H

#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "nonConformal/pointPatches/nonConformalCyclic/nonConformalCyclicPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalCyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalCyclicPointPatchField
:
    public pointPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(nonConformalCyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalCyclic field
        nonConformalCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, nonConformalCyclic field and dictionary
        nonConformalCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalCyclicPointPatchField
        (
            const nonConformalCyclicPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting nonConformalCyclic field reference
        nonConformalCyclicPointPatchField
        (
            const nonConformalCyclicPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalCyclicPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting nonConformalCyclic field
        //  reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalCyclicPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformal/pointPatchFields/nonConformalCyclic/nonConformalCyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
