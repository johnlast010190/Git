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
    (c) 2017 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "probes/probes.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note:chould check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::probes::sampleAndWrite(const VolField<Type>& vField)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        if (newFileFormat_)
        {
            Info<< "Writing in new format" << endl;
            columnatedFileWriter* fileWriter = fileWriters_[vField.name()];
            fileWriter->writeTime();
            forAll(values, probei)
            {
                if (coorFramePtr_)
                {
                    vector probeLocation = this->operator[](probei);
                    if (localOnOutput_)
                    {
                        probeLocation =
                            coorFramePtr_->coorSys().localPosition
                            (
                                probeLocation
                            );
                    }
                    fileWriter->writeDelimited(probeLocation);
                }
                fileWriter->writeDelimited(values[probei]);
            }
            fileWriter->endLine();
        }
        else
        {
            unsigned int w = IOstream::defaultPrecision() + 7;
            OFstream& os = *probeFilePtrs_[vField.name()];

            os  << setw(w) << Time::timeName(vField.time().timeOutputValue());

            forAll(values, probei)
            {
                if (includeOutOfBounds_ || processor_[probei] != -1)
                {
                    if (coorFramePtr_)
                    {
                        vector probeLocation = this->operator[](probei);
                        if (localOnOutput_)
                        {
                            probeLocation =
                                coorFramePtr_->coorSys().localPosition
                                (
                                    probeLocation
                                );
                        }
                        os  << ' ' << setw(w) << probeLocation;
                    }
                    os  << ' ' << setw(w) << values[probei];
                }
            }
            os  << endl;
        }
    }
}


template<class Type>
void Foam::probes::sampleAndWrite(const SurfaceField<Type>& sField)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        if (newFileFormat_)
        {
            columnatedFileWriter* fileWriter = fileWriters_[sField.name()];
            fileWriter->writeTime();
            forAll(values, probei)
            {
                if (coorFramePtr_)
                {
                    vector probeLocation = this->operator[](probei);
                    if (localOnOutput_)
                    {
                        probeLocation =
                            coorFramePtr_->coorSys().localPosition
                            (
                                probeLocation
                            );
                    }
                    fileWriter->writeDelimited(probeLocation);
                }
                fileWriter->writeDelimited(values[probei]);
            }
            fileWriter->endLine();
        }
        else
        {
            unsigned int w = IOstream::defaultPrecision() + 7;
            OFstream& os = *probeFilePtrs_[sField.name()];

            os  << setw(w) << Time::timeName(sField.time().timeOutputValue());

            forAll(values, probei)
            {
                if (includeOutOfBounds_ || processor_[probei] != -1)
                {
                    if (coorFramePtr_)
                    {
                        vector probeLocation = this->operator[](probei);
                        if (localOnOutput_)
                        {
                            probeLocation =
                                coorFramePtr_->coorSys().localPosition
                                (
                                    probeLocation
                                );
                        }
                        os  << ' ' << setw(w) << probeLocation;
                    }
                    os  << ' ' << setw(w) << values[probei];
                }
            }
            os  << endl;
        }
    }
}


template<class Type>
void Foam::probes::sampleAndWrite(const fieldGroup<Type>& fields)
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
                        obr(),
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
            objectRegistry::const_iterator iter = obr().find(fields[fieldi]);

            if (iter.found() && iter()->type() == VolField<Type>::typeName)
            {
                sampleAndWrite(lookupObject<VolField<Type>>(fields[fieldi]));
            }
        }
    }
}


template<class Type>
void Foam::probes::sampleAndWriteSurfaceFields(const fieldGroup<Type>& fields)
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
                        obr(),
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
            objectRegistry::const_iterator iter = obr().find(fields[fieldi]);

            if
            (
                iter.found()
             && iter()->type()
             == SurfaceField<Type>::typeName
            )
            {
                sampleAndWrite
                (
                    lookupObject<SurfaceField<Type>>(fields[fieldi])
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const VolField<Type>& vField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);
    tmp<Field<Type>> tValues(new Field<Type>(this->size(), unsetVal));
    Field<Type>& values = tValues.ref();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type>> interpolator
        (
            interpolation<Type>::New(interpolationScheme_, vField)
        );

        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                const vector& position = operator[](probei);

                values[probei] = interpolator().interpolate
                (
                    position,
                    elementList_[probei],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                values[probei] = vField[elementList_[probei]];
            }
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::probes::sample(const word& fieldName) const
{
    return sample(lookupObject<VolField<Type>>(fieldName));
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::probes::sample
(
    const SurfaceField<Type>& sField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);
    tmp<Field<Type>> tValues(new Field<Type>(this->size(), unsetVal));
    Field<Type>& values = tValues.ref();

    forAll(*this, probei)
    {
        if (faceList_[probei] >= 0)
        {
            values[probei] = sField[faceList_[probei]];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sampleSurfaceFields(const word& fieldName) const
{
    return sample(lookupObject<SurfaceField<Type>>(fieldName));
}

// ************************************************************************* //
