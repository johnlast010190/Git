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
    (c) 2015 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "parFvFieldReconstructor.H"
#include "db/Time/Time.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "fields/fvsPatchFields/constraint/empty/emptyFvsPatchField.H"
#include "db/IOobjectList/IOobjectList.H"
#include "meshes/polyMesh/polyDistributionMap/polyDistributionMap.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"

#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"
#include "distributedUnallocatedDirectFieldMapper.H"
#include "distributedUnallocatedDirectFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const DimensionedField<Type, volMesh>& fld
) const
{
    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<Type> internalField(mapper(fld));

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<DimensionedField<Type, volMesh>> tfield
    (
        new DimensionedField<Type, volMesh>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField
        )
    );

    tfield.ref().oriented() = fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    DimensionedField<Type, volMesh> fld
    (
        fieldIoObject,
        procMesh_
    );

    // Distribute onto baseMesh
    return reconstructFvVolumeInternalField(fld);
}


// Reconstruct a field onto the baseMesh
template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::parFvFieldReconstructor::reconstructFvVolumeField
(
    const VolField<Type>& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<Type> internalField(mapper(fld.internalField()));


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<fvPatchField<Type>> patchFields(fld.mesh().boundary().size());

    const typename VolField<Type>::Boundary&
        bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        if (patchFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<fvPatchField<Type>> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const fvPatchField<Type>& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    pfld,
                    basePatch,
                    DimensionedField<Type, volMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(basePatchFields, patchI)
    {
        if (patchI >= patchFields.size() || !patchFields.set(patchI))
        {
            basePatchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    emptyFvPatchField<Type>::typeName,
                    baseMesh_.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
    }

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<VolField<Type>> tfield
    (
        new VolField<Type>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );

    tfield.ref().oriented()= fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::parFvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    VolField<Type> fld
    (
        fieldIoObject,
        procMesh_,
        true,
        false //central boundaries
    );

    // Distribute onto baseMesh
    return reconstructFvVolumeField(fld);
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::parFvFieldReconstructor::reconstructFvSurfaceField
(
    const SurfaceField<Type>& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap()
    );

    // Create flat field of internalField + all patch fields
    Field<Type> flatFld(fld.mesh().nFaces(), Type(Zero));
    SubList<Type>(flatFld, fld.internalField().size()) = fld.internalField();
    forAll(fld.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& fvp = fld.boundaryField()[patchI];

        SubList<Type>(flatFld, fvp.size(), fvp.patch().start()) = fvp;
    }

    // Map all faces
    Field<Type> internalField(mapper(flatFld, fld.oriented()()));

    // Trim to internal faces (note: could also have special mapper)
    internalField.setSize
    (
        min
        (
            internalField.size(),
            baseMesh_.nInternalFaces()
        )
    );


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<fvsPatchField<Type>> patchFields(fld.mesh().boundary().size());

    const typename SurfaceField<Type>::Boundary&
        bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        if (patchFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<fvsPatchField<Type>> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const fvsPatchField<Type>& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    pfld,
                    basePatch,
                    DimensionedField<Type, surfaceMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(basePatchFields, patchI)
    {
        if (patchI >= patchFields.size() || !patchFields.set(patchI))
        {
            basePatchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    emptyFvsPatchField<Type>::typeName,
                    baseMesh_.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
    }

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<SurfaceField<Type>> tfield
    (
        new SurfaceField<Type>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );

    tfield.ref().oriented() = fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::parFvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    SurfaceField<Type> fld
    (
        fieldIoObject,
        procMesh_,
        true,
        false //central boundaries
    );

    return reconstructFvSurfaceField(fld);
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvVolumeInternalFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName = DimensionedField<Type, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                tmp<DimensionedField<Type, volMesh>> tfld
                (
                    reconstructFvVolumeInternalField<Type>(*fieldIter())
                );

                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvVolumeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        VolField<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            const word& name = fieldIter()->name();

            if
            (
                (selectedFields.empty() || selectedFields.found(name))
             && name != "cellDist"
            )
            {
                Info<< "        " << name << endl;

                tmp<VolField<Type>> tfld
                (
                    reconstructFvVolumeField<Type>(*fieldIter())
                );
                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvSurfaceFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        SurfaceField<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                tmp<SurfaceField<Type>> tfld
                (
                    reconstructFvSurfaceField<Type>(*fieldIter())
                );
                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


// ************************************************************************* //
