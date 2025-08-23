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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
void Foam::meshDistribute::printFieldInfo(const fvMesh& mesh)
{
    const UPtrList<GeoField> fields(mesh.fields<GeoField>());

    forAll(fields, i)
    {
        const GeoField& field = fields[i];

        Pout<< "Field:" << field.name() << " internal size:" << field.size()
            << endl;

        forAll(field.boundaryField(), patchi)
        {
            Pout<< "    " << patchi
                << ' ' << field.boundaryField()[patchi].patch().name()
                << ' ' << field.boundaryField()[patchi].type()
                << ' ' << field.boundaryField()[patchi].size()
                << endl;
        }
    }
}


template<class T, class Mesh>
void Foam::meshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, T>>& bfields,
    const fvMesh& mesh
) const
{
    // Save whole boundary field

    typedef GeometricField<T, fvsPatchField, Mesh> fieldType;

    const UPtrList<fieldType> fields(mesh.fields<fieldType>());

    bfields.setSize(fields.size());

    forAll(fields, i)
    {
        const fieldType& field = fields[i];

        bfields.set(i, field.boundaryField().clone().ptr());
    }
}


template<class T, class Mesh>
void Foam::meshDistribute::mapBoundaryFields
(
    const polyTopoChangeMap& map,
    const PtrList<FieldField<fvsPatchField, T>>& oldBfields,
    fvMesh& mesh
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, Mesh> fieldType;

    UPtrList<fieldType> fields(mesh.fields<fieldType>());

    forAll(fields, i)
    {
        fieldType& field = fields[i];
        typename fieldType::Boundary& bfield = field.boundaryFieldRef();

        const FieldField<fvsPatchField, T>& oldBfield = oldBfields[i];

        // Pull from old boundary field into bfield

        forAll(bfield, patchi)
        {
            fvsPatchField<T>& patchField = bfield[patchi];
            label facei = patchField.patch().start();

            forAll(patchField, i)
            {
                label oldFacei = faceMap[facei++];

                // Find patch and local patch face oldFacei was in.
                forAll(oldPatchStarts, oldPatchi)
                {
                    label oldLocalI = oldFacei - oldPatchStarts[oldPatchi];

                    if
                    (
                        oldLocalI >= 0
                     && oldLocalI < oldBfield[oldPatchi].size()
                    )
                    {
                        patchField[i] = oldBfield[oldPatchi][oldLocalI];
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::meshDistribute::saveInternalFields(PtrList<Field<T>>& ifields) const
{
    const UPtrList<SurfaceField<T>> fields(mesh_.fields<SurfaceField<T>>());

    ifields.setSize(fields.size());

    forAll(fields, i)
    {
        const SurfaceField<T>& field = fields[i];

        ifields.set(i, field.primitiveField().clone());
    }
}


template<class T>
void Foam::meshDistribute::mapExposedFaces
(
    const polyTopoChangeMap& map,
    const PtrList<Field<T>>& oldFields
)
{
    // Set boundary values of exposed internal faces

    const labelList& faceMap = map.faceMap();

    UPtrList<SurfaceField<T>> fields(mesh_.fields<SurfaceField<T>>());

    forAll(fields, i)
    {
        SurfaceField<T>& field = fields[i];
        //const bool oriented = field.oriented()();

        typename SurfaceField<T>::Boundary& bfield = field.boundaryFieldRef();
        const bool negateIfFlipped = isFlux(field);

        const Field<T>& oldInternal = oldFields[i];

        // Pull from old internal field into bfield

        forAll(bfield, patchi)
        {
            fvsPatchField<T>& patchField = bfield[patchi];

            forAll(patchField, i)
            {
                const label faceI = patchField.patch().start() + i;
                const label oldFaceI = faceMap[faceI];

                if (oldFaceI < oldInternal.size())
                {
                    patchField[i] = oldInternal[oldFaceI];

                    if (negateIfFlipped && map.flipFaceFlux().found(faceI))
                    {
                        patchField[i] = flipOp()(patchField[i]);
                    }
                }
            }
        }
    }
}


template<class GeoField, class PatchFieldType>
void Foam::meshDistribute::initPatchFields
(
    const typename GeoField::value_type& initVal,
    fvMesh& mesh
)
{
    // Init patch fields of certain type

    UPtrList<GeoField> fields(mesh.fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        forAll(bfield, patchi)
        {
            if (isA<PatchFieldType>(bfield[patchi]))
            {
                bfield[patchi].forceAssign(initVal);
            }
        }
    }
}


template<class GeoField>
void Foam::meshDistribute::correctBoundaryConditions()
{
    // CorrectBoundaryConditions patch fields of certain type

    UPtrList<GeoField> fields(mesh_.fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];
        field.correctBoundaryConditions();
    }
}


template<class GeoField>
void Foam::meshDistribute::resetUpdate()
{
    UPtrList<GeoField> fields(mesh_.fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                field.boundaryFieldRef()[patchi].resetUpdate();
            }
        }
    }
}


template<class GeoField>
void Foam::meshDistribute::sendFields
(
    const label domain,
    const wordList& fieldNames,
    const fvMeshSubset& subsetter,
    Ostream& toNbr
)
{
    // Send fields. Note order supplied so we can receive in exactly the same
    // order.
    // Note that field gets written as entry in dictionary so we
    // can construct from subdictionary.
    // (since otherwise the reading as-a-dictionary mixes up entries from
    // consecutive fields)
    // The dictionary constructed is:
    //  volScalarField
    //  {
    //      p {internalField ..; boundaryField ..;}
    //      k {internalField ..; boundaryField ..;}
    //  }
    //  volVectorField
    //  {
    //      U {internalField ...  }
    //  }

    // volVectorField {U {internalField ..; boundaryField ..;}}

    toNbr << GeoField::typeName << token::NL << token::BEGIN_BLOCK << token::NL;
    forAll(fieldNames, i)
    {
        if (debug)
        {
            Pout<< "Subsetting field " << fieldNames[i]
                << " for domain:" << domain << endl;
        }

        // Send all fieldNames. This has to be exactly the same set as is
        // being received!
        const GeoField& field =
            subsetter.baseMesh().lookupObject<GeoField>(fieldNames[i]);

        // Suppress warning "mapper does not map all values."
        tmp<GeoField> tsubfield =
            subsetter.interpolate(field, false); // debug true issue!

        toNbr
            << fieldNames[i] << token::NL << token::BEGIN_BLOCK
            << tsubfield
            << token::NL << token::END_BLOCK << token::NL;
    }
    toNbr << token::END_BLOCK << token::NL;
}


template<class GeoField>
void Foam::meshDistribute::initFields
(
    fvMesh& baseMesh,
    fvMesh& mesh,
    const autoPtr<fvMeshSubset>& subsetterPtr,
    PtrList<GeoField>& iFields
)
{
    const UPtrList<GeoField> fields(baseMesh.fields<GeoField>());

    iFields.setSize(fields.size());

    forAll(fields, i)
    {
        const GeoField& field = fields[i];

        const word& name = field.name();

        // Create zero sized field and send
        if (subsetterPtr.valid())
        {
            if (debug)
            {
                Info<< "Initializing zero sized field: " << name
                    << endl;
            }

            tmp<GeoField> tsubfield = subsetterPtr().interpolate(field);

            OStringStream ostr(IOstream::BINARY);
            ostr<< tsubfield();

            IStringStream istr
            (
                ostr.str(),
                IOstream::BINARY
            );

            dictionary fieldDict(istr);

            iFields.set
            (
                i,
                new GeoField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    fieldDict
                )
            );
        }
    }
}


template<class GeoField>
void Foam::meshDistribute::receiveFields
(
    const label domain,
    const wordList& fieldNames,
    typename GeoField::Mesh& mesh,
    PtrList<GeoField>& fields,
    const dictionary& fieldDicts
)
{
    // Opposite of sendFields

    if (debug)
    {
        Pout<< "Receiving fields " << fieldNames
            << " from domain:" << domain << endl;
    }

    fields.setSize(fieldNames.size());

    forAll(fieldNames, i)
    {
        if (debug)
        {
            Pout<< "Constructing field " << fieldNames[i]
                << " from domain:" << domain << endl;
        }

        fields.set
        (
            i,
            new GeoField
            (
                IOobject
                (
                    fieldNames[i],
                    mesh.thisDb().time().timeName(),
                    mesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fieldDicts.subDict(fieldNames[i])
            )
        );
    }
}


// ************************************************************************* //
