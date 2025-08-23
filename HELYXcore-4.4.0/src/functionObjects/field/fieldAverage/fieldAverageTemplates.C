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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldAverage/fieldAverage.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "surfFields/surfFields/surfFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAverage::readMeanFieldType(const label fieldi)
{
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& baseFieldName = faItems_[fieldi].fieldName();

    IOobject meanFieldIO
    (
        meanFieldName,
        obr().time().timeName(),
        obr(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if
    (
        meanFieldIO.headerOk()
     && meanFieldIO.headerClassName() == Type::typeName
    )
    {
        if (obr().found(meanFieldName))
        {
            Log << "    Cannot read average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging for field." << endl;

            faItems_[fieldi].mean() = false;
        }
        else
        {
            Log << "    Reading field " << meanFieldName << endl;

            const Type& baseField = lookupObject<Type>(baseFieldName);

            obr().store(new Type(meanFieldIO, baseField));
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::readMeanField(const label fieldi)
{
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (faItems_[fieldi].mean())
    {
        readMeanFieldType<VolField<Type>>(fieldi);
        readMeanFieldType<SurfaceField<Type>>(fieldi);
        readMeanFieldType<SurfFieldType>(fieldi);
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::initialiseMeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    if (obr().foundObject<Type>(meanFieldName))
    {
        // Do nothing
    }
    else if (obr().found(meanFieldName))
    {
        Log << "    Cannot initialise average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].mean() = false;
    }
    else
    {
        Log << "    Initialising field " << meanFieldName << endl;

        const Type& baseField = obr().lookupObject<Type>(fieldName);

        obr().store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr().time().timeName(),
                    obr()
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::initialiseMeanField
(
    const label fieldi
)
{
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (faItems_[fieldi].mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (obr().foundObject<VolField<Type>>(fieldName))
        {
            initialiseMeanFieldType<VolField<Type>>(fieldi);
        }
        else if (obr().foundObject<SurfaceField<Type>>(fieldName))
        {
            initialiseMeanFieldType<SurfaceField<Type>>(fieldi);
        }
        else if (obr().foundObject<SurfFieldType>(fieldName))
        {
            initialiseMeanFieldType<SurfFieldType>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::readPrime2MeanFieldType
(
    const label fieldi
)
{
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();
    const word& baseFieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    IOobject prime2MeanFieldIO
    (
        prime2MeanFieldName,
        obr().time().timeName(),
        obr(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if
    (
        prime2MeanFieldIO.headerOk()
     && prime2MeanFieldIO.headerClassName() == Type2::typeName
    )
    {
        if (obr().found(prime2MeanFieldName))
        {
            Log << "    Cannot read average field " << prime2MeanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging for field." << endl;

            faItems_[fieldi].mean() = false;
        }
        else
        {
            Log << "    Reading field " << prime2MeanFieldName << endl;

            const Type1& baseField = lookupObject<Type1>(baseFieldName);
            const Type1& meanField = lookupObject<Type1>(meanFieldName);

            obr().store
            (
                new Type2(prime2MeanFieldIO, sqr(baseField) - sqr(meanField))
            );
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::readPrime2MeanField
(
    const label fieldi
)
{
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << faItems_[fieldi].fieldName() << nl
                << exit(FatalError);
        }

        readPrime2MeanFieldType
        <VolField<Type1>, VolField<Type2>>(fieldi);

        readPrime2MeanFieldType
        <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);

        readPrime2MeanFieldType
        <SurfFieldType1, SurfFieldType2>(fieldi);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    if (obr().foundObject<Type2>(prime2MeanFieldName))
    {
        // Do nothing
    }
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot initialise average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        Log << "    Initialising field " << prime2MeanFieldName << nl;

        const Type1& baseField = obr().lookupObject<Type1>(fieldName);
        const Type1& meanField = obr().lookupObject<Type1>(meanFieldName);

        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(),
                    obr()
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanField
(
    const label fieldi
)
{
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr().foundObject<VolField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldType
            <VolField<Type1>, VolField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfaceField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldType
            <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfFieldType1>(fieldName))
        {
            initialisePrime2MeanFieldType
            <SurfFieldType1, SurfFieldType2>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldWelType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    if (obr().foundObject<Type2>(prime2MeanFieldName))
    {
        // Do nothing
    }
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot initialise average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        Log << "    Initialising field " << prime2MeanFieldName << nl;

        const Type1& baseField = obr().lookupObject<Type1>(fieldName);
        const Type1& meanField = obr().lookupObject<Type1>(meanFieldName);

        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(),
                    obr()
                ),
                (baseField - meanField)*(baseField - meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldWel
(
    const label fieldi
)
{
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr().foundObject<VolField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldWelType
            <VolField<Type1>, VolField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfaceField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldWelType
            <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfFieldType1>(fieldName))
        {
            initialisePrime2MeanFieldWelType
            <SurfFieldType1, SurfFieldType2>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldWelSymmType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    if (obr().foundObject<Type2>(prime2MeanFieldName))
    {
        // Do nothing
    }
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot initialise average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        Log << "    Initialising field " << prime2MeanFieldName << nl;

        const Type1& baseField = obr().lookupObject<Type1>(fieldName);
        const Type1& meanField = obr().lookupObject<Type1>(meanFieldName);

        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(),
                    obr()
                ),
                sqr(baseField - meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::initialisePrime2MeanFieldWelSymm
(
    const label fieldi
)
{
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr().foundObject<VolField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldWelSymmType
            <VolField<Type1>, VolField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfaceField<Type1>>(fieldName))
        {
            initialisePrime2MeanFieldWelSymmType
            <SurfaceField<Type1>, SurfaceField<Type2>>(fieldi);
        }
        else if (obr().foundObject<SurfFieldType1>(fieldName))
        {
            initialisePrime2MeanFieldWelSymmType
            <SurfFieldType1, SurfFieldType2>(fieldi);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::copyFieldType
(
    const word& fieldName,
    const word& outputFieldName,
    bool& activationFlag,
    bool& operationFlag
)
{
    // Field has been found, so set active flag to true
    activationFlag = true;

    Log << "    Reading/initialising field " << outputFieldName << endl;

    if (foundObject<Type>(outputFieldName))
    {}
    else if (obr().found(outputFieldName))
    {
        Log << "    Cannot allocate average field " << outputFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        operationFlag = false;
    }
    else
    {
        const Type& baseField = lookupObject<Type>(fieldName);

        // Store on registry
        obr().store
        (
            new Type
            (
                IOobject
                (
                    outputFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::addMinField(const label fieldi)
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (faItems_[fieldi].minFlag())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (foundObject<VolFieldType>(fieldName))
        {
            copyFieldType<VolFieldType>
            (
                fieldName,
                faItems_[fieldi].minFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].minFlag()
            );
        }
        else if (foundObject<SurfaceFieldType>(fieldName))
        {
            copyFieldType<SurfaceFieldType>
            (
                fieldName,
                faItems_[fieldi].minFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].minFlag()
            );
        }
        else if (foundObject<SurfFieldType>(fieldName))
        {
            copyFieldType<SurfFieldType>
            (
                fieldName,
                faItems_[fieldi].minFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].minFlag()
            );
        }
    }
}

template<class Type>
void Foam::functionObjects::fieldAverage::addMaxField(const label fieldi)
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (faItems_[fieldi].maxFlag())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (foundObject<VolFieldType>(fieldName))
        {
            copyFieldType<VolFieldType>
            (
                fieldName,
                faItems_[fieldi].maxFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].maxFlag()
            );
        }
        else if (foundObject<SurfaceFieldType>(fieldName))
        {
            copyFieldType<SurfaceFieldType>
            (
                fieldName,
                faItems_[fieldi].maxFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].maxFlag()
            );
        }
        else if (foundObject<SurfFieldType>(fieldName))
        {
            copyFieldType<SurfFieldType>
            (
                fieldName,
                faItems_[fieldi].maxFieldName(),
                faItems_[fieldi].active(),
                faItems_[fieldi].maxFlag()
            );
        }
    }
}

template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    const Type& baseField = lookupObject<Type>(fieldName);

    Type& meanField =
        lookupObjectRef<Type>(faItems_[fieldi].meanFieldName());

    scalar dt = obr().time().deltaTValue();
    scalar Dt = totalTime_[fieldi];

    if (faItems_[fieldi].iterBase())
    {
        dt = 1;
        Dt = scalar(totalIter_[fieldi]);
    }

    scalar beta = dt/Dt;

    if (faItems_[fieldi].window() > 0)
    {
        const scalar w = faItems_[fieldi].window();

        if (Dt - dt >= w)
        {
            beta = dt/w;
        }
    }

    if (faItems_[fieldi].timeIntegral() && faItems_[fieldi].timeBase())
    {
        if (totalIter_[fieldi] == 1)
        {
            meanField = baseField*dt;
        }
        else
        {
            meanField = meanField + baseField*dt;
        }
    }
    else
    {
        meanField = (1 - beta)*meanField + beta*baseField;
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].fieldName();

            if (obr().foundObject<VolField<Type>>(fieldName))
            {
                calculateMeanFieldType<VolField<Type>>(fieldi);
            }
            else if (obr().foundObject<SurfaceField<Type>>(fieldName))
            {
                calculateMeanFieldType<SurfaceField<Type>>(fieldi);
            }
            else if (obr().foundObject<SurfFieldType>(fieldName))
            {
                calculateMeanFieldType<SurfFieldType>(fieldi);
            }
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        const Type1& meanField =
            lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

        scalar dt = obr().time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (faItems_[fieldi].iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;

        if (faItems_[fieldi].window() > 0)
        {
            const scalar w = faItems_[fieldi].window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        prime2MeanField =
            (1 - beta)*prime2MeanField
          + beta*sqr(baseField)
          - sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculateMeanPrime2MeanFieldWelType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        Type1& meanField =
            lookupObjectRef<Type1>(faItems_[fieldi].meanFieldName());

        scalar dt = obr().time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (faItems_[fieldi].iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;
        if (faItems_[fieldi].window() > 0)
        {
            const scalar w = faItems_[fieldi].window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        Type1 meanFieldOld = meanField;

        if (faItems_[fieldi].timeIntegral() && faItems_[fieldi].timeBase())
        {
            if (totalIter_[fieldi] == 1)
            {
                meanField = baseField*dt;
            }
            else
            {
                meanField = meanFieldOld + baseField*dt;
            }
        }
        else
        {
            meanField = (1 - beta)*meanFieldOld + beta*baseField;
        }

        if (faItems_[fieldi].prime2Mean())
        {
            Type2& prime2MeanField =
                lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());
            prime2MeanField = (1 - beta)*prime2MeanField
                + beta*((baseField - meanFieldOld)*(baseField - meanField));
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::
calculateMeanPrime2MeanFieldWelSymmType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        Type1& meanField =
            lookupObjectRef<Type1>(faItems_[fieldi].meanFieldName());

        scalar dt = obr().time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (faItems_[fieldi].iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;
        if (faItems_[fieldi].window() > 0)
        {
            const scalar w = faItems_[fieldi].window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        Type1 meanFieldOld = meanField;

        if (faItems_[fieldi].timeIntegral() && faItems_[fieldi].timeBase())
        {
            if (totalIter_[fieldi] == 1)
            {
                meanField = baseField*dt;
            }
            else
            {
                meanField = meanFieldOld + baseField*dt;
            }
        }
        else
        {
            meanField = (1 - beta)*meanFieldOld + beta*baseField;
        }

        if (faItems_[fieldi].prime2Mean())
        {
            Type2& prime2MeanField =
                lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());
            prime2MeanField = (1 - beta)*prime2MeanField
                + beta*symm
                    (
                        ((baseField - meanFieldOld)*(baseField - meanField))
                    );
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFields() const
{
    typedef VolField<Type1> VolFieldType1;
    typedef SurfaceField<Type1> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef VolField<Type2> VolFieldType2;
    typedef SurfaceField<Type2> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            calculatePrime2MeanFieldType<VolFieldType1, VolFieldType2>(fieldi);
            calculatePrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                fieldi
            );
            calculatePrime2MeanFieldType<SurfFieldType1, SurfFieldType2>
            (
                fieldi
            );
        }
    }
}

template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculateMeanPrime2MeanFieldsWel()
const
{
    typedef VolField<Type1> VolFieldType1;
    typedef SurfaceField<Type1> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef VolField<Type2> VolFieldType2;
    typedef SurfaceField<Type2> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            calculateMeanPrime2MeanFieldWelType
                <VolFieldType1, VolFieldType2>(fieldi);
            calculateMeanPrime2MeanFieldWelType
                <SurfaceFieldType1, SurfaceFieldType2>(fieldi);
            calculateMeanPrime2MeanFieldWelType<SurfFieldType1, SurfFieldType2>
            (
                fieldi
            );
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculateMeanPrime2MeanFieldsWelSymm()
const
{
    typedef VolField<Type1> VolFieldType1;
    typedef SurfaceField<Type1> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef VolField<Type2> VolFieldType2;
    typedef SurfaceField<Type2> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            calculateMeanPrime2MeanFieldWelSymmType
                <VolFieldType1, VolFieldType2>(fieldi);
            calculateMeanPrime2MeanFieldWelSymmType
                <SurfaceFieldType1, SurfaceFieldType2>(fieldi);
            calculateMeanPrime2MeanFieldWelSymmType
                <SurfFieldType1, SurfFieldType2>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2MeanType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& meanField =
            lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

        prime2MeanField += sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef VolField<Type1> VolFieldType1;
    typedef SurfaceField<Type1> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef VolField<Type2> VolFieldType2;
    typedef SurfaceField<Type2> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].prime2Mean())
        {
            addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>(fieldi);
            addMeanSqrToPrime2MeanType<SurfaceFieldType1, SurfaceFieldType2>
            (
                fieldi
            );
            addMeanSqrToPrime2MeanType<SurfFieldType1, SurfFieldType2>(fieldi);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateExtFieldType
(
    const word& fieldName,
    const word& outputFieldName,
    scalar nx
) const
{
    if (foundObject<Type>(fieldName))
    {
        const Type& baseField = lookupObject<Type>(fieldName);

        Type& outputField = const_cast<Type&>
        (
            lookupObject<Type>(outputFieldName)
        );

        outputField = nx*max(nx*outputField, nx*baseField);
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMinFields() const
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].minFlag())
        {
            calculateExtFieldType<VolFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].minFieldName(),
                -1.0
            );
            calculateExtFieldType<SurfaceFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].minFieldName(),
                -1.0
            );
            calculateExtFieldType<SurfFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].minFieldName(),
                -1.0
            );
        }
    }
}

template<class Type>
void Foam::functionObjects::fieldAverage::calculateMaxFields() const
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].maxFlag())
        {
            calculateExtFieldType<VolFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].maxFieldName(),
                1.0
            );
            calculateExtFieldType<SurfaceFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].maxFieldName(),
                1.0
            );
            calculateExtFieldType<SurfFieldType>
            (
                faItems_[fieldi].fieldName(),
                faItems_[fieldi].maxFieldName(),
                1.0
            );
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFieldType
(
    const word& fieldName
) const
{
    if (foundObject<Type>(fieldName))
    {
        const Type& f = lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFields() const
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& fieldName = faItems_[fieldi].meanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
        if (faItems_[fieldi].prime2Mean())
        {
            const word& fieldName = faItems_[fieldi].prime2MeanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
        if (faItems_[fieldi].minFlag())
        {
            const word& fieldName = faItems_[fieldi].minFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
        if (faItems_[fieldi].maxFlag())
        {
            const word& fieldName = faItems_[fieldi].maxFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
    }
}


// ************************************************************************* //
