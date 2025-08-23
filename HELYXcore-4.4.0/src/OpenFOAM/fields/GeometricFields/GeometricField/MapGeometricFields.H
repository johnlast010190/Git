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
    (c) 2025 Engys Ltd.

Class
    Foam::MapInternalField

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MapGeometricFields_H
#define MapGeometricFields_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper, class GeoMesh>
class MapInternalField
{
public:

    MapInternalField()
    {}

    void operator()
    (
        DimensionedField<Type, GeoMesh>& field,
        const MeshMapper& mapper
    ) const;
};


//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations
//  for mapping of internal fields depending on mesh type.
template
<
    class Type,
    template<class> class PatchField,
    class MeshMapper,
    class GeoMesh
>
void MapGeometricFields
(
    const MeshMapper& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

    HashTable<const FieldType*> fields
    (
        mapper.thisDb().objectRegistry::template
            lookupClass<FieldType>()
    );

    forAllConstIter(typename HashTable<const FieldType*>, fields, fieldIter)
    {
        if (GeoMesh::Mesh::geometryFields.found(fieldIter()->name())) continue;

        FieldType& field = const_cast<FieldType&>(*fieldIter());

        if (&field.mesh() == &mapper.mesh())
        {
            if (polyMesh::debug)
            {
                Info<< "Mapping " << field.typeName << ' ' << field.name()
                    << endl;
            }

            // Map the internal field
            MapInternalField<Type, MeshMapper, GeoMesh>()(field.ref(), mapper);

            // Map the patch fields
            typename FieldType::Boundary& bfield = field.boundaryFieldRef();

            forAll(bfield, patchi)
            {
                // Cannot check sizes for patch fields because of
                // empty fields in FV and because point fields get their size
                // from the patch which has already been resized
                //

                bfield[patchi].autoMap(mapper.boundaryMap()[patchi]);
            }

            field.instance() = field.time().timeName();
        }
        else if (polyMesh::debug)
        {
            Info<< "Not mapping " << field.typeName << ' ' << field.name()
                << " since originating mesh differs from that of mapper."
                << endl;
        }
    }
}


template<class GeoField>
void AddPatchFields
(
    objectRegistry& obr,
    const label insertPatchi,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const typename GeoField::value_type& defaultPatchValue
)
{
    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());

    const wordList fldNames(flds.sortedToc());

    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        if (bfld.size() != fld.mesh().boundary().size())
        {
            FatalErrorInFunction << "bfld size:" << bfld.size()
                << " mesh size:" << fld.mesh().boundary().size()
                << exit(FatalError);
        }

        if (patchFieldDict.found(fld.name()))
        {
            bfld.set
            (
                insertPatchi,
                GeoField::Patch::New
                (
                    fld.mesh().boundary()[insertPatchi],
                    fld(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
        else
        {
            bfld.set
            (
                insertPatchi,
                GeoField::Patch::New
                (
                    defaultPatchFieldType,
                    fld.mesh().boundary()[insertPatchi],
                    fld()
                )
            );
            bfld[insertPatchi].forceAssign(defaultPatchValue);
        }
    }
}


template<class GeoField>
void ReorderPatchFields
(
    objectRegistry& obr,
    const labelUList& newToOld
)
{
    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());
    const wordList fldNames(flds.sortedToc());
    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];
        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        bfld.shuffle(newToOld);
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void UpdateGIBFields
(
    polyMesh& mesh
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {

        typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

        if (GeoMesh::Mesh::geometryFields.found(fieldIter()->name())) continue;

        FieldType& field = const_cast<FieldType&>(*fieldIter());

        // Map the patch fields
        typename FieldType::Boundary& bfield = field.boundaryFieldRef();

        forAll(bfield, patchi)
        {
            // Cannot check sizes for patch fields because of
            // empty fields in FV and because point fields get their size
            // from the patch which has already been resized
            //
            bfield[patchi].updateGIB();
        }
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void StoreGIBFields
(
    polyMesh& mesh
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {

        typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

        if (GeoMesh::Mesh::geometryFields.found(fieldIter()->name())) continue;

        FieldType& field = const_cast<FieldType&>(*fieldIter());

        // Map the patch fields
        typename FieldType::Boundary& bfield = field.boundaryFieldRef();

        forAll(bfield, patchi)
        {
            bfield[patchi].storeGIB();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
