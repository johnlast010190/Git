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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::nonConformalDiscreteMixingPointPatchField

Description
    Constraint pointPatchField for non-conformal discrete mixing couplings.

SourceFiles
    nonConformalDiscreteMixingPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingPointPatchField_H
#define nonConformalDiscreteMixingPointPatchField_H

#include "fields/pointPatchFields/basic/coupled/coupledPointPatchField.H"
#include "nonConformal/pointPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPointPatch.H"
#include "interpolations/primitivePatchInterpolation/PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nonConformalDiscreteMixingPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalDiscreteMixingPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private Data

        //- Local reference cast into the non-conformal point patch
        const nonConformalDiscreteMixingPointPatch& nonConformalPatch_;

        //- Owner side patch interpolation pointer
        mutable autoPtr<PrimitivePatchInterpolation<primitivePatch>> ppiPtr_;

        //- Neighbour side patch interpolation pointer
        mutable autoPtr<PrimitivePatchInterpolation<primitivePatch>> nbrPpiPtr_;


    // Private Member Functions

        //- Return owner side patch interpolation
        const PrimitivePatchInterpolation<primitivePatch>& ppi() const
        {
            if (!ppiPtr_.valid())
            {
                ppiPtr_.reset
                (
                    new PrimitivePatchInterpolation<primitivePatch>
                    (
                        nonConformalPatch_.nonConformalPolyPatch()
                    )
                );
            }

            return ppiPtr_();
        }

        //- Return neighbour side patch interpolation
        const PrimitivePatchInterpolation<primitivePatch>& nbrPpi() const
        {
            if (!nbrPpiPtr_.valid())
            {
                nbrPpiPtr_.reset
                (
                    new PrimitivePatchInterpolation<primitivePatch>
                    (
                        nonConformalPatch_.nonConformalPolyPatch().nbrPatch()
                    )
                );
            }

            return nbrPpiPtr_();
        }


public:

    //- Runtime type information
    TypeName(nonConformalDiscreteMixingPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalDiscreteMixingPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalDiscreteMixingPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalDiscreteMixingPointPatchField
        (
            const nonConformalDiscreteMixingPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalDiscreteMixingPointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonConformalDiscreteMixingPointPatchField
        (
            const nonConformalDiscreteMixingPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone, setting internal field reference.
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalDiscreteMixingPointPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        ) override
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformal/pointPatchFields/nonConformalDiscreteMixing/nonConformalDiscreteMixingPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
