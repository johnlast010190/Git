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

#include "fields/pointPatchFields/derived/fixedNormalSlip/fixedNormalSlipPointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    slipPointPatchField<Type>(p, iF),
    n_(vector::max)
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    slipPointPatchField<Type>(p, iF, dict),
    n_(dict.lookup("n"))
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const fixedNormalSlipPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    slipPointPatchField<Type>(ptf, p, iF, mapper),
    n_(ptf.n_)
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const fixedNormalSlipPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    slipPointPatchField<Type>(ptf, iF),
    n_(ptf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedNormalSlipPointPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    tmp<Field<Type>> tvalues =
        transform(I - n_*n_, this->patchInternalField());

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInternalField(iF, tvalues());
}


template<class Type>
void Foam::fixedNormalSlipPointPatchField<Type>::write(Ostream& os) const
{
    slipPointPatchField<Type>::write(os);
    os.writeEntry("n", n_);
}


// ************************************************************************* //
