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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/GeometricField/GeometricIOBoundaryField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const BoundaryMesh& bmesh
)
:
    regIOobject(io),
    GeoBoundaryField(bmesh)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const BoundaryMesh& bmesh,
    const Internal& field,
    const word& patchFieldType,
    const wordList& actualPatchTypes
)
:
    regIOobject(io),
    GeoBoundaryField(bmesh, field, patchFieldType, actualPatchTypes)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const BoundaryMesh& bmesh,
    const Internal& field,
    const wordList& patchFieldTypes,
    const wordList& constraintTypes
)
:
    regIOobject(io),
    GeoBoundaryField(bmesh, field, patchFieldTypes, constraintTypes)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const BoundaryMesh& bmesh,
    const Internal& field,
    const PtrList<PatchField<Type>>& ptfl
)
:
    regIOobject(io),
    GeoBoundaryField(bmesh, field, ptfl)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const Internal& field,
    const GeoBoundaryField& btf
)
:
    regIOobject(io),
    GeoBoundaryField(field, btf)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const Internal& field,
    GeoBoundaryField& btf,
    bool reuse
)
:
    regIOobject(io),
    GeoBoundaryField(field, btf, reuse)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const GeometricIOBoundaryField<Type, PatchField, GeoMesh>& btf
)
:
    regIOobject(btf),
    GeoBoundaryField(btf)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricIOBoundaryField<Type, PatchField, GeoMesh>::GeometricIOBoundaryField
(
    const IOobject& io,
    const BoundaryMesh& bmesh,
    const Internal& field,
    const dictionary& dict
)
:
    regIOobject(io),
    GeoBoundaryField(bmesh, field, dict)
{}


// ************************************************************************* //
