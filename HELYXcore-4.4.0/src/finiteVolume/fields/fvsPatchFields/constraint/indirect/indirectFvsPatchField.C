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
    (c) 2011 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvsPatchFields/constraint/indirect/indirectFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
indirectFvsPatchField<Type>::indirectFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF)
{}


template<class Type>
indirectFvsPatchField<Type>::indirectFvsPatchField
(
    const indirectFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
indirectFvsPatchField<Type>::indirectFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
indirectFvsPatchField<Type>::indirectFvsPatchField
(
    const indirectFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>(ptf)
{}


template<class Type>
indirectFvsPatchField<Type>::indirectFvsPatchField
(
    const indirectFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class Type>
tmp<Field<Type>> indirectFvsPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );
}


template<class Type>
tmp<Field<Type>> indirectFvsPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    //return *this;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type>> indirectFvsPatchField<Type>::gradientInternalCoeffs() const
{
    //return -pTraits<Type>::one*this->patch().deltaCoeffs();
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type>> indirectFvsPatchField<Type>::gradientBoundaryCoeffs() const
{
    //return this->patch().deltaCoeffs()*(*this);
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}
*/

template<class Type>
tmp<Field<Type>> indirectFvsPatchField<Type>::gradientBoundaryValue() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
void indirectFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
