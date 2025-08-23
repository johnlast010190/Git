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
    (c) 2015-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents() const
{
    return pow
    (
        this->solutionD(),
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents2() const
{
    return pow
    (
        (this->solutionD() + Vector<label>::one)/2,
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class Type>
void Foam::fvMesh::stabiliseEmptyDirections(Field<Type>& t) const
{
    const typename pTraits<Type>::labelType v = validComponents2<Type>();
    const Type ident = pTraits<Type>::I;
    for (label cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        if (!v[cmpt])
        {
            t.replace(cmpt, ident[cmpt]);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fvMesh::stabiliseEmptyDirections
(
    GeometricField<Type, PatchField, GeoMesh>& t
) const
{
    const typename pTraits<Type>::labelType v = validComponents2<Type>();
    for (label cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        if (!v[cmpt])
        {
            const dimensioned<typename pTraits<Type>::cmptType> dimCmpt
            (
                "I", t.dimensions(), pTraits<Type>::I[cmpt]
            );
            t.replace(cmpt, dimCmpt);
        }
    }
}


template<class Type>
Foam::labelHashSet Foam::fvMesh::findActivePatchIDs() const
{
    const polyBoundaryMesh& bm = this->boundaryMesh();
    labelHashSet patchIDs(bm.findPatchIDs<Type>());

    labelHashSet patchIDsMod(patchIDs.size());

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        label patchi =  iter.key();
        if (this->boundary()[patchi].isActive())
        {
            patchIDsMod.insert(patchi);
        }
    }

    return patchIDs;
}


template<class GeoField>
Foam::UPtrList<GeoField> Foam::fvMesh::fields(const bool strict) const
{
    HashTable<GeoField*> fields
    (
        const_cast<fvMesh&>(*this).lookupClass<GeoField>(strict)
    );
    UPtrList<GeoField> curFields(fields.size());

    label i = 0;
    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (!geometryFields.found(iter()->name()))
        {
            curFields.set(i++, iter());
        }
    }
    curFields.setSize(i);

    return curFields;
}


template<class GeoField>
Foam::UPtrList<GeoField> Foam::fvMesh::curFields() const
{
    HashTable<GeoField*> fields
    (
        const_cast<fvMesh&>(*this).lookupClass<GeoField>()
    );
    UPtrList<GeoField> curFields(fields.size());

    label i = 0;
    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (!geometryFields.found(iter()->name()) && !iter()->isOldTime())
        {
            curFields.set(i++, iter());
        }
    }
    curFields.setSize(i);

    return curFields;
}


template<class Type, template<class> class GeoField>
void Foam::fvMesh::storeOldTimeFields()
{
    UPtrList<GeoField<Type>> curFields(this->curFields<GeoField<Type>>());

    forAll(curFields, i)
    {
        curFields[i].storeOldTimes();
    }
}


template<template<class> class GeoField>
void Foam::fvMesh::storeOldTimeFields()
{
    #define StoreOldTimeFields(Type, nullArg) \
        storeOldTimeFields<Type, GeoField>();

    FOR_ALL_FIELD_TYPES(StoreOldTimeFields);

    #undef StoreOldTimeFields
}


template<class Type, template<class> class GeoField>
void Foam::fvMesh::nullOldestTimeFields()
{
    UPtrList<GeoField<Type>> curFields(this->curFields<GeoField<Type>>());

    forAll(curFields, i)
    {
        curFields[i].nullOldestTime();
    }
}


template<template<class> class GeoField>
void Foam::fvMesh::nullOldestTimeFields()
{
    #define nullOldestTimeFields(Type, nullArg) \
        nullOldestTimeFields<Type, GeoField>();

    FOR_ALL_FIELD_TYPES(nullOldestTimeFields);

    #undef nullOldestTimeFields
}


// ************************************************************************* //
