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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/cyclic/cyclicFvPatchField.H"
#include "fields/Fields/transformField/transformField.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFvPatchField<Type>::interpPatchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicPatch().nbrPatch().faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = transform().transform(iField[nbrFaceCells[facei]]);
    }

    return tpnf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isA<cyclicFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    this->evaluate(Pstream::commsTypes::blocking);
}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isA<cyclicFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (this->calculatedMode())
    {
        // The boundary values were stored as the patch neighbour field values
        // in evaluate(). We use those stored values instead of the
        // implementation of patchNeighbourField below because a derived type
        // may have overridden it, but will be substituted by this constraint
        // type when used in field calculations
        return *this;
    }
    else
    {
        return interpPatchNeighbourField();
    }
}


template<class Type>
const Foam::cyclicFvPatchField<Type>&
Foam::cyclicFvPatchField<Type>::nbrPatchField() const
{
    const VolField<Type>& fld =
    static_cast<const VolField<Type>&>
    (
        this->primitiveField()
    );

    return refCast<const cyclicFvPatchField<Type>>
    (
        fld.boundaryField()[this->cyclicPatch().nbrPatchID()]
    );
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // We co-opt the patch field to store the neighbour field to be used
    // in field calculations. We can do this since it isn't used for snGrad
    // (which uses the coupled value), and other operators that use the face
    // value will use the coupled interpolated value, rather than the
    // patch value.
    if (this->calculatedMode())
    {
        // This is an exception: If evaluate is called on this base class in
        // calculated mode (which shouldn't really happen, but is allowed
        // for convenience), we avoid a circular assignment
        Field<Type>::operator=(interpPatchNeighbourField());
    }
    else
    {
        Field<Type>::operator=
        (
            this->patchNeighbourField(commsType)
        );
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicPatch().nbrPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicPatch().nbrPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    NotImplemented;
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
}


// ************************************************************************* //
