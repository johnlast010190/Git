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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshDistribute::printFieldInfo(const fvMesh& mesh)
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
void Foam::fvMeshDistribute::saveBoundaryFields
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
void Foam::fvMeshDistribute::mapBoundaryFields
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
void Foam::fvMeshDistribute::initMapExposedFaces
(
    PtrList<Field<T>>& ifields
) const
{
    const UPtrList<SurfaceField<T>> fields(mesh_.fields<SurfaceField<T>>());

    ifields.setSize(fields.size());

    forAll(fields, i)
    {
        ifields.set(i, Field<T>(mesh_.nFaces()));

        const SurfaceField<T>& field = fields[i];

        SubList<T>(ifields[i], field.primitiveField().size()) =
            field.primitiveField();

        forAll(field.boundaryField(), patchi)
        {
            const fvsPatchField<T>& pfield = field.boundaryField()[patchi];

            SubList<T>(ifields[i], pfield.size(), pfield.patch().start()) =
                pfield;
        }
    }
}


template<class T>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const polyTopoChangeMap& map,
    const PtrList<Field<T>>& oldFields
)
{
    UPtrList<SurfaceField<T>> fields(mesh_.fields<SurfaceField<T>>());

    forAll(fields, i)
    {
        SurfaceField<T>& field = fields[i];

        const Field<T>& oldField = oldFields[i];

        const bool negateIfFlipped = isFlux(field);

        forAll(field.boundaryField(), patchi)
        {
            fvsPatchField<T>& patchField = field.boundaryFieldRef()[patchi];

            forAll(patchField, i)
            {
                const label facei = patchField.patch().start() + i;
                const label oldFacei = map.faceMap()[facei];

                if (oldFacei < map.nOldInternalFaces())
                {
                    if (negateIfFlipped && map.flipFaceFlux().found(facei))
                    {
                        patchField[i] = flipOp()(oldField[oldFacei]);
                    }
                    else
                    {
                        patchField[i] = oldField[oldFacei];
                    }
                }
                else
                {
                    patchField[i] = oldField[oldFacei];
                }
            }
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::correctCoupledPatchFields(fvMesh& mesh)
{
    UPtrList<GeoField> fields(mesh.fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(bfield, patchi)
            {
                if (bfield[patchi].coupled())
                {
                    bfield[patchi].initEvaluate(Pstream::defaultCommsType);
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(bfield, patchi)
            {
                if (bfield[patchi].coupled())
                {
                    bfield[patchi].evaluate(Pstream::defaultCommsType);
                }
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule =
                mesh_.globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                if (bfield[patchEvali].coupled())
                {
                    if (patchSchedule[patchEvali].init)
                    {
                        bfield[patchSchedule[patchEvali].patch]
                            .initEvaluate(Pstream::commsTypes::scheduled);
                    }
                    else
                    {
                        bfield[patchSchedule[patchEvali].patch]
                            .evaluate(Pstream::commsTypes::scheduled);
                    }
                }
            }
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::resetUpdate()
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
void Foam::fvMeshDistribute::sendFields
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
void Foam::fvMeshDistribute::initFields
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
void Foam::fvMeshDistribute::receiveFields
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
