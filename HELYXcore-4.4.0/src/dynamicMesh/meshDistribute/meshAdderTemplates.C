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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"
#include "fvMeshAdder/directPointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::meshAdder::MapVolField
(
    const mapAddedPolyMesh& meshMap,

    VolField<Type>& fld,
    const VolField<Type>& fldToAdd,
    const bool fullyMapped
)
{
    const fvMesh& mesh = fld.mesh();

    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        Field<Type> oldInternalField(fld.primitiveField());

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nCells());

        intFld.rmap(oldInternalField, meshMap.oldCellMap());
        intFld.rmap(fldToAdd.primitiveField(), meshMap.addedCellMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    typename VolField<Type>::
    Boundary& bfld = fld.boundaryFieldRef();

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchStarts = meshMap.oldPatchStarts();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1              // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Override mapping (for use in e.g. fvMeshDistribute where
                // it sorts mapping out itself)
                if (fullyMapped)
                {
                    patchMapper.hasUnmapped() = false;
                }

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - fld.boundaryField()[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchi)
        {
            label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1          // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    // Override mapping (for use in e.g. fvMeshDistribute where
                    // it sorts mapping out itself)
                    if (fullyMapped)
                    {
                        patchMapper.hasUnmapped() = false;
                    }

                    bfld.set
                    (
                        newPatchi,
                        fvPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi], // added field
                            mesh.boundary()[newPatchi],       // new fvPatch
                            fld(),   // new int. field
                            patchMapper                       // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFacei = oldPatch.start()+i;
                        label newFacei = meshMap.addedFaceMap()[addedFacei];
                        label patchFacei = newFacei-newPatch.start();
                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::meshAdder::MapVolFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd,
    const bool fullyMapped
)
{
    HashTable<const VolField<Type>*> fields
    (
        mesh.objectRegistry::lookupClass
        <VolField<Type>>
        ()
    );

    HashTable<const VolField<Type>*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass
        <VolField<Type>>
        ()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        DebugPout
            << "MapVolFields : Storing old time for " << fieldIter()->name()
            << endl;

        const_cast<VolField<Type>*>(fieldIter())
            ->storeOldTimes();
    }


    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (fvMesh::geometryFields.found(fieldIter()->name())) continue;

        VolField<Type>& fld =
            const_cast<VolField<Type>&>
            (
                *fieldIter()
            );

        if (fieldsToAdd.found(fld.name()))
        {
            const VolField<Type>& fldToAdd =
                *fieldsToAdd[fld.name()];

            DebugPout
                << "MapVolFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapVolField<Type>(meshMap, fld, fldToAdd, fullyMapped);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::meshAdder::MapSurfaceField
(
    const mapAddedPolyMesh& meshMap,

    SurfaceField<Type>& fld,
    const SurfaceField<Type>& fldToAdd,
    const bool fullyMapped
)
{
    const fvMesh& mesh = fld.mesh();
    const labelList& oldPatchStarts = meshMap.oldPatchStarts();

    typename SurfaceField<Type>::
    Boundary& bfld = fld.boundaryFieldRef();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nInternalFaces());

        intFld.rmap(oldField, meshMap.oldFaceMap());
        intFld.rmap(fldToAdd, meshMap.addedFaceMap());


        // Faces that were boundary faces but are not anymore.
        // Use owner value (so lowest numbered cell, i.e. from 'old' not 'added'
        // mesh)
        forAll(bfld, patchi)
        {
            const fvsPatchField<Type>& pf = bfld[patchi];

            label start = oldPatchStarts[patchi];

            forAll(pf, i)
            {
                label newFacei = meshMap.oldFaceMap()[start + i];

                if (newFacei >= 0 && newFacei < mesh.nInternalFaces())
                {
                    intFld[newFacei] = pf[i];
                }
            }
        }
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1      // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Override mapping (for use in e.g. fvMeshDistribute where
                // it sorts mapping out itself)
                if (fullyMapped)
                {
                    patchMapper.hasUnmapped() = false;
                }

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - bfld[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvsPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchi)
        {
            label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1                  // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    // Override mapping (for use in e.g. fvMeshDistribute where
                    // it sorts mapping out itself)
                    if (fullyMapped)
                    {
                        patchMapper.hasUnmapped() = false;
                    }

                    bfld.set
                    (
                        newPatchi,
                        fvsPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi],// added field
                            mesh.boundary()[newPatchi],      // new fvPatch
                            fld(),  // new int. field
                            patchMapper                      // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFacei = oldPatch.start()+i;
                        label newFacei = meshMap.addedFaceMap()[addedFacei];
                        label patchFacei = newFacei-newPatch.start();
                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::meshAdder::MapSurfaceFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd,
    const bool fullyMapped
)
{
    typedef SurfaceField<Type> fldType;

    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>()
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        DebugPout
            << "MapSurfaceFields : Storing old time for "
            << fieldIter()->name() << endl;

        const_cast<fldType*>(fieldIter())->storeOldTimes();
    }


    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (fvMesh::geometryFields.found(fieldIter()->name())) continue;

        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            DebugPout
                << "MapSurfaceFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapSurfaceField<Type>(meshMap, fld, fldToAdd, fullyMapped);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::meshAdder::MapPointField
(
    const pointMesh& mesh,
    const mapAddedPolyMesh& meshMap,
    const labelListList& oldMeshPoints,

    GeometricField<Type, pointPatchField, pointMesh>& fld,
    const GeometricField<Type, pointPatchField, pointMesh>& fldToAdd
)
{
    // This is a bit tricky:
    // - mesh pointed to by fld is invalid
    // - pointPatches pointed to be fld are invalid

    typename GeometricField<Type, pointPatchField, pointMesh>::
    Boundary& bfld = fld.boundaryFieldRef();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.size());

        intFld.rmap(oldField, meshMap.oldPointMap());
        intFld.rmap(fldToAdd, meshMap.addedPointMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundary().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                const labelList& oldMp = oldMeshPoints[patchi];
                const pointPatch& newPp = mesh.boundary()[newPatchi];
                const labelList& newMeshPoints = newPp.meshPoints();
                Map<label> newMeshPointMap(2*newMeshPoints.size());
                forAll(newMeshPoints, ppi)
                {
                    newMeshPointMap.insert(newMeshPoints[ppi], ppi);
                }

                labelList newToOld(newPp.size(), -1);
                forAll(oldMp, oldPointi)
                {
                    label newPointi = oldMp[oldPointi];

                    Map<label>::const_iterator fnd =
                        newMeshPointMap.find(newPointi);
                    if (fnd == newMeshPointMap.end())
                    {
                        // Possibly an internal point
                    }
                    else
                    {
                        // Part of new patch
                        newToOld[fnd()] = oldPointi;
                    }
                }

                directPointPatchFieldMapper patchMapper(newToOld);

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - bfld[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    pointPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new pointPatch
                        fld(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchi)
        {
            label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const pointPatch& oldPatch = fldToAdd.mesh().boundary()[patchi];
                const labelList& oldMp = oldPatch.meshPoints();
                const pointPatch& newPp = mesh.boundary()[newPatchi];
                const labelList& newMeshPoints = newPp.meshPoints();
                Map<label> newMpm(2*newMeshPoints.size());
                forAll(newMeshPoints, ppi)
                {
                    newMpm.insert(newMeshPoints[ppi], ppi);
                }

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    labelList newToAdded(newPp.size(), -1);
                    forAll(oldMp, oldPointi)
                    {
                        label newPointi = oldMp[oldPointi];
                        Map<label>::const_iterator fnd = newMpm.find(newPointi);
                        if (fnd == newMpm.end())
                        {
                            // Possibly an internal point
                        }
                        else
                        {
                            // Part of new patch
                            newToAdded[fnd()] = oldPointi;
                        }
                    }

                    directPointPatchFieldMapper patchMapper(newToAdded);

                    bfld.set
                    (
                        newPatchi,
                        pointPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi],// added field
                            mesh.boundary()[newPatchi],      // new pointPatch
                            fld(),  // new int. field
                            patchMapper                      // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList oldToNew(oldPatch.size(), -1);
                    forAll(oldMp, oldPointi)
                    {
                        label newPointi = oldMp[oldPointi];
                        Map<label>::const_iterator fnd = newMpm.find(newPointi);
                        if (fnd != newMpm.end())
                        {
                            // Part of new patch
                            oldToNew[oldPointi] = fnd();
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        oldToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::meshAdder::MapPointFields
(
    const mapAddedPolyMesh& meshMap,
    const pointMesh& mesh,
    const labelListList& oldMeshPoints,
    const objectRegistry& meshToAdd
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fldType;

    HashTable<const fldType*> fields(mesh.thisDb().lookupClass<fldType>());
    HashTable<const fldType*> fieldsToAdd(meshToAdd.lookupClass<fldType>());

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (debug)
        {
            Pout<< "MapPointFields : Storing old time for "
                << fieldIter()->name() << endl;
        }

        const_cast<fldType*>(fieldIter())->storeOldTimes();
    }


    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (fvMesh::geometryFields.found(fieldIter()->name())) continue;

        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            if (debug)
            {
                Pout<< "MapPointFields : mapping " << fld.name()
                    << " and " << fldToAdd.name() << endl;
            }

            MapPointField<Type>(mesh, meshMap, oldMeshPoints, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::meshAdder::MapDimField
(
    const mapAddedPolyMesh& meshMap,

    DimensionedField<Type, volMesh>& fld,
    const DimensionedField<Type, volMesh>& fldToAdd
)
{
    const fvMesh& mesh = fld.mesh();

    // Store old field
    Field<Type> oldField(fld);

    fld.setSize(mesh.nCells());

    fld.rmap(oldField, meshMap.oldCellMap());
    fld.rmap(fldToAdd, meshMap.addedCellMap());
}


template<class Type>
void Foam::meshAdder::MapDimFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
)
{
    typedef DimensionedField<Type, volMesh> fldType;

    // Note: use strict flag on lookupClass to avoid picking up
    //       volFields
    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>(true)
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>(true)
    );

    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            DebugPout
                << "MapDimFields : mapping " << fld.name()
                << " and " << fldToAdd.name() << endl;

            MapDimField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


// ************************************************************************* //
