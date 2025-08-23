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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2010-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "db/dictionary/dictionary.H"
#include "fvMesh/fvMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/distributedFvPatchFieldMapper.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    Field<Type>(p.size(), pTraits<Type>::zero),
    patch_(p),
    internalField_(iF),
    patchType_(word::null),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    patchType_(word::null),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null)),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_
    (
        dict.lookupOrDefault<Type>
        (
            "defaultGIBValue",
            pTraits<Type>::zero
        )
    )
{
    if (valueRequired)
    {
        if (dict.found("value"))
        {
            fvsPatchField<Type>::operator=
            (
                Field<Type>("value", dict, p.size())
            );
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "essential 'value' entry not provided"
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{
    if (mappingRequired)
    {
        mapper(*this, ptf);
    }
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField(const fvsPatchField<Type>& ptf)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(ptf.defaultGIBValue_)
{}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(ptf.defaultGIBValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fvsPatchField<Type>::db() const
{
    if (notNull(internalField_))
    {
        return internalField_.db();
    }
    else
    {
        return patch_.boundaryMesh().mesh();
    }
}


template<class Type>
void Foam::fvsPatchField<Type>::check(const fvsPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorInFunction
            << "different patches for fvsPatchField<Type>s"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvsPatchField<Type>::autoMap(const fvPatchFieldMapper& m)
{
    // Pass oriented flag to fvPatchFieldMappers supporting
    // distributed functionality
    if (isA<distributedFvPatchFieldMapper>(m))
    {
        const bool oriented = internalField_.oriented()();

        const distributedFvPatchFieldMapper& dm =
            dynamic_cast<const distributedFvPatchFieldMapper&>(m);

        dm(*this, *this, oriented);
    }
    else
    {
        m(*this, *this);
    }
}


template<class Type>
void Foam::fvsPatchField<Type>::rmap
(
    const fvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvsPatchField<Type>::reset(const fvsPatchField<Type>& ptf)
{
    Field<Type>::reset(ptf);
}

template<class Type>
void Foam::fvsPatchField<Type>::updateGIB()
{
    const polyPatch& pthis = this->patch().patch();
    if (isA<indirectPolyPatch>(pthis))
    {
        storeGIB();
        const indirectPolyPatch& dpp =
            refCast<const indirectPolyPatch>(pthis);
        this->setSize(dpp.size());
        fvsPatchField<Type>::forceAssign(pTraits<Type>::zero);
    }
}

template<class Type>
void Foam::fvsPatchField<Type>::autoMapGIB(const gibFvPatchFieldMapper&)
{}

template<class Type>
void Foam::fvsPatchField<Type>::storeGIB()
{
    const polyPatch& pthis = this->patch().patch();
    if (isA<indirectPolyPatch>(pthis))
    {
        oldTimeFieldPtr_.reset(new Field<Type>(*this));
    }
}


template<class Type>
void Foam::fvsPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (patchType_.size())
    {
        os.writeEntry("patchType", patchType_);
    }

    // The 'value' field is only written by patch fields that need it
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const fvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const fvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const fvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::fvsPatchField<Type>::forceAssign
(
    const fvsPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::forceAssign
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::forceAssign
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvsPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFieldNew.C"

// ************************************************************************* //
