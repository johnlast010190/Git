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

\*---------------------------------------------------------------------------*/

#include "probes/patchProbes.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const VolField<Type>& vField
)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& probeStream = *probeFilePtrs_[vField.name()];

        probeStream
            << setw(w)
            << Time::timeName(vField.time().timeOutputValue());

        forAll(values, probei)
        {
            probeStream << ' ' << setw(w) << values[probei];
        }
        probeStream << endl;
    }
}


template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const SurfaceField<Type>& sField
)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& probeStream = *probeFilePtrs_[sField.name()];

        probeStream
            << setw(w)
            << Time::timeName(sField.time().timeOutputValue());

        forAll(values, probei)
        {
            probeStream << ' ' << setw(w) << values[probei];
        }
        probeStream << endl;
    }
}


template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const fieldGroup<Type>& fields
)
{
    forAll(fields, fieldi)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                VolField<Type>
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fields[fieldi]);

            if
            (
                iter.found()
             && iter()->type()
             == VolField<Type>::typeName
            )
            {
                sampleAndWrite
                (
                    mesh_.lookupObject
                    <VolField<Type>>
                    (
                        fields[fieldi]
                    )
                );
            }
        }
    }
}


template<class Type>
void Foam::patchProbes::sampleAndWriteSurfaceFields
(
    const fieldGroup<Type>& fields
)
{
    forAll(fields, fieldi)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                SurfaceField<Type>
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fields[fieldi]);

            if
            (
                iter.found()
             && iter()->type()
             == SurfaceField<Type>::typeName
            )
            {
                sampleAndWrite
                (
                    mesh_.lookupObject
                    <SurfaceField<Type>>
                    (
                        fields[fieldi]
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample
(
    const VolField<Type>& vField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    tmp<Field<Type>> tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = vField.boundaryField()[patchi][localFacei];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<VolField<Type>>
        (
            fieldName
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample
(
    const SurfaceField<Type>& sField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    tmp<Field<Type>> tValues
    (
        new Field<Type>(this->size(), unsetVal)
    );

    Field<Type>& values = tValues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = sField.boundaryField()[patchi][localFacei];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


// ************************************************************************* //
