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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    referenceFrameFvPatch<Type>(p, iF)
{}


template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchField<Type>(p, iF, value),
    referenceFrameFvPatch<Type>(p, iF)
{}


template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<Type>(p, iF, dict, valueRequired),
    referenceFrameFvPatch<Type>(dict, p, iF)
{}


template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    referenceFrameFvPatch<Type>
    (
        p,
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        mapper,
        iF
    )
{}


template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fixedValueFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        ptf.internalField()
    )
{}


template<class Type>
Foam::fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
    const fixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        iF
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    referenceFrameFvPatch<Type>::autoMap(m);
}


template<class Type>
void Foam::fixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);
    const fixedValueFvPatchField<Type>& fvvptf =
        refCast<const fixedValueFvPatchField<Type>>(ptf);
    referenceFrameFvPatch<Type>::rmap(fvvptf, addr);
}


template<class Type>
void Foam::fixedValueFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchField<Type>::autoMapGIB(mapper);
    referenceFrameFvPatch<Type>::autoMapGIB(mapper);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>(new Field<Type>(this->size(), Zero));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void Foam::fixedValueFvPatchField<Type>::frameFieldUpdate
(
    Field<Type>& tf,
    bool setInValue,
    bool setToGlobal
)
{
    this->addFrameVelocity(tf, setInValue, setToGlobal);
    this->forceAssign(tf);
}


template<class Type>
void Foam::fixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    referenceFrameFvPatch<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
