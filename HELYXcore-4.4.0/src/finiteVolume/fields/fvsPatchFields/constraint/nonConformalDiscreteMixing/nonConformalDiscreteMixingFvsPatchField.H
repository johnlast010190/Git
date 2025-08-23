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
    Foam::nonConformalDiscreteMixingFvsPatchField

Description
    Foam::nonConformalDiscreteMixingFvsPatchField

SourceFiles
    nonConformalDiscreteMixingFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingFvsPatchField_H
#define nonConformalDiscreteMixingFvsPatchField_H

#include "fields/fvsPatchFields/basic/coupled/coupledFvsPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonConformalDiscreteMixingFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalDiscreteMixingFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the non-conformal discrete-mixing patch
        const nonConformalDiscreteMixingFvPatch& nonConformalPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalDiscreteMixingFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalDiscreteMixingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalDiscreteMixingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalDiscreteMixingFvsPatchField
        //  onto a new patch
        nonConformalDiscreteMixingFvsPatchField
        (
            const nonConformalDiscreteMixingFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nonConformalDiscreteMixingFvsPatchField
        (
            const nonConformalDiscreteMixingFvsPatchField<Type>&
        );

        //- Copy constructor setting field reference
        nonConformalDiscreteMixingFvsPatchField
        (
            const nonConformalDiscreteMixingFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const override
        {
            return tmp<fvsPatchField<Type>>
            (
                new nonConformalDiscreteMixingFvsPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const override
        {
            return tmp<fvsPatchField<Type>>
            (
                new nonConformalDiscreteMixingFvsPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes commsType
            ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvsPatchFields/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
