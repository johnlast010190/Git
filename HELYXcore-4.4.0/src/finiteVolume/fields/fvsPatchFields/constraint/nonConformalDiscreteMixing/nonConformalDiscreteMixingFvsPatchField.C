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

#include "fields/fvsPatchFields/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::
nonConformalDiscreteMixingFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingFvPatch>(p))
{}


template<class Type>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::
nonConformalDiscreteMixingFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingFvPatch>(p))
{
    if (!isA<nonConformalDiscreteMixingFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::
nonConformalDiscreteMixingFvsPatchField
(
    const nonConformalDiscreteMixingFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingFvPatch>(p))
{
    if (!isA<nonConformalDiscreteMixingFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::
nonConformalDiscreteMixingFvsPatchField
(
    const nonConformalDiscreteMixingFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    nonConformalPatch_(ptf.nonConformalPatch_)
{}


template<class Type>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::
nonConformalDiscreteMixingFvsPatchField
(
    const nonConformalDiscreteMixingFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    nonConformalPatch_(ptf.nonConformalPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::nonConformalDiscreteMixingFvsPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    typedef SurfaceField<Type> geoField;
    const geoField& gf = refCast<const geoField>(this->internalField());

    const nonConformalDiscreteMixingFvPatch& ncdmp =
        refCast<const nonConformalDiscreteMixingFvPatch>(this->patch());

    const Field<Type>& nbf = gf.boundaryField()[ncdmp.nbrPatchID()];

    tmp<Field<Type>> tpnf = nonConformalPatch_.interpolate(nbf);

    ncdmp.transform().transform(tpnf.ref(), tpnf());

    return tpnf;
}


// ************************************************************************* //
