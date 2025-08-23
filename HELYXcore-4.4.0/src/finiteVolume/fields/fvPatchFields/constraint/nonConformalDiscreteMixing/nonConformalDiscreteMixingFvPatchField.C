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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::
nonConformalDiscreteMixingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    nonConformalDiscreteMixingPatch_
    (
        refCast<const nonConformalDiscreteMixingFvPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::
nonConformalDiscreteMixingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    nonConformalDiscreteMixingPatch_
    (
        refCast<const nonConformalDiscreteMixingFvPatch>(p)
    )
{
    this->evaluate(Pstream::commsTypes::blocking);
}


template<class Type>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::
nonConformalDiscreteMixingFvPatchField
(
    const nonConformalDiscreteMixingFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    nonConformalDiscreteMixingPatch_
    (
        refCast<const nonConformalDiscreteMixingFvPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::
nonConformalDiscreteMixingFvPatchField
(
    const nonConformalDiscreteMixingFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    nonConformalDiscreteMixingPatch_(ptf.nonConformalDiscreteMixingPatch_)
{}


template<class Type>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::
nonConformalDiscreteMixingFvPatchField
(
    const nonConformalDiscreteMixingFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    nonConformalDiscreteMixingPatch_(ptf.nonConformalDiscreteMixingPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::nonConformalDiscreteMixingFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        nonConformalPatch().nbrPatch().faceCells();

    Field<Type> nif(iField, nbrFaceCells);

    tmp<Field<Type>> tpnf = nonConformalDiscreteMixingPatch_.interpolate(nif);

    return tpnf;
}


template<class Type>
void Foam::nonConformalDiscreteMixingFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=(this->patchNeighbourField(commsType));

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::nonConformalDiscreteMixingFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells = nonConformalPatch().nbrPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    pnf = nonConformalPatch().interpolate(pnf);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::nonConformalDiscreteMixingFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells = nonConformalPatch().nbrPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    pnf = nonConformalPatch().interpolate(pnf);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::nonConformalDiscreteMixingFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
}


// ************************************************************************* //
