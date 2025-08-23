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
    (c) 2022-2023 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Description
    Perform mapping of finite volume fields required by stitching.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/fvMeshStitchers/fvMeshStitcher/conformedFvsPatchField.H"
#include "nonConformal/boundary/nonConformalBoundary.H"
#include "fvMesh/fvPatches/constraint/nonConformal/nonConformalFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    UPtrList<GeoField<Type>> fields(mesh_.fields<GeoField<Type>>());

    forAll(fields, i)
    {
        forAll(mesh_.boundary(), patchi)
        {
            typename GeoField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (isA<nonConformalFvPatch>(pf.patch()))
            {
                pf.autoMap(fvPatchFieldSetSizer(pf.patch().size()));
            }
        }
    }
}


template<template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    #define ResizePatchFields(Type, nullArg) \
        resizePatchFields<Type, GeoField>();

    FOR_ALL_FIELD_TYPES(ResizePatchFields);

    #undef ResizePatchFields

    //- Resize vectorN-type fields
    resizePatchFields<vector4, GeoField>();
}


template<class Type>
void Foam::fvMeshStitcher::preConformSurfaceFields()
{
    UPtrList<SurfaceField<Type>> fields(mesh_.curFields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];

        for (label ti = 0; ti <= field.nOldTimes(false); ti++)
        {
            conformedFvsPatchField<Type>::conform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformSurfaceFields()
{
    UPtrList<SurfaceField<Type>> fields(mesh_.curFields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];

        for (label ti = 0; ti <= field.nOldTimes(false); ti++)
        {
            conformedFvsPatchField<Type>::unconform
            (
                boundaryFieldRefNoUpdate(field.oldTime(ti))
            );
        }
    }
}


template<class Type>
void Foam::fvMeshStitcher::evaluateVolFields()
{
    UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const label nReq = Pstream::nRequests();

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (isA<nonConformalFvPatch>(pf.patch()))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                boundaryFieldRefNoUpdate(fields[i])[patchi];

            if (isA<nonConformalFvPatch>(pf.patch()))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class GeoField>
typename GeoField::Boundary& Foam::fvMeshStitcher::boundaryFieldRefNoUpdate
(
    GeoField& fld
)
{
    return const_cast<typename GeoField::Boundary&>(fld.boundaryField());
}


// ************************************************************************* //
