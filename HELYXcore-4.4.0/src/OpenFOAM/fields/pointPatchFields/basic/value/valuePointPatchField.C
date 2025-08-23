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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2012 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/pointPatchFields/basic/value/valuePointPatchField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::valuePointPatchField<Type>::checkFieldSize() const
{
    if (this->size() != this->patch().size())
    {
        FatalErrorInFunction
            << "field does not correspond to patch. " << endl
            << "Field size: " << size() << " patch size: "
            << this->patch().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::valuePointPatchField<Type>::valuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(p, iF),
    Field<Type>(p.size())
{}


template<class Type>
Foam::valuePointPatchField<Type>::valuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    pointPatchField<Type>(p, iF, dict),
    Field<Type>(p.size())
{
    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else if (!valueRequired)
    {
        Field<Type>::operator=(Zero);
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::valuePointPatchField<Type>::valuePointPatchField
(
    const valuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    pointPatchField<Type>(ptf, p, iF, mapper),
    Field<Type>(mapper(ptf))
{}


template<class Type>
Foam::valuePointPatchField<Type>::valuePointPatchField
(
    const valuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(ptf, iF),
    Field<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::valuePointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    m(*this, *this);
}


template<class Type>
void Foam::valuePointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap
    (
        refCast<const valuePointPatchField<Type>>
        (
            ptf
        ),
        addr
    );
}


template<class Type>
void Foam::valuePointPatchField<Type>::reset
(
    const pointPatchField<Type>& ptf
)
{
    Field<Type>::reset(refCast<const valuePointPatchField<Type>>(ptf));
}


template<class Type>
void Foam::valuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInternalField(iF, *this);

    pointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::valuePointPatchField<Type>::updateFromInternalField()
{
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    tmp<Field<Type>> valuesPtr = this->patchInternalField(iF);
    Field<Type> &values = valuesPtr.ref();

    this->transfer(values);

    pointPatchField<Type>::updateFromInternalField();
}


template<class Type>
void Foam::valuePointPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInternalField(iF, *this);

    pointPatchField<Type>::evaluate();
}


template<class Type>
void Foam::valuePointPatchField<Type>::write(Ostream& os) const
{
    pointPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::valuePointPatchField<Type>::operator=
(
    const valuePointPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::valuePointPatchField<Type>::operator=
(
    const pointPatchField<Type>& ptf
)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void Foam::valuePointPatchField<Type>::operator=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::valuePointPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::valuePointPatchField<Type>::forceAssign
(
    const valuePointPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::valuePointPatchField<Type>::forceAssign
(
    const pointPatchField<Type>& ptf
)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void Foam::valuePointPatchField<Type>::forceAssign
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::valuePointPatchField<Type>::forceAssign
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// ************************************************************************* //
