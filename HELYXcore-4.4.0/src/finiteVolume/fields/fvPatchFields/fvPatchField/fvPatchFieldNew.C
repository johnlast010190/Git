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
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << " : " << p.type()
            << endl;
    }

    tmp<fvPatchField<Type>> tfvp;

    const auto ctor =
        ctorTableLookup
        (
            "patchField type", patchConstructorTable_(), patchFieldType
        );
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTable_().find(p.type());

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (patchTypeCstrIter != patchConstructorTable_().end())
        {
            tfvp = patchTypeCstrIter->second(p, iF);
        }
        else
        {
            tfvp = ctor(p, iF);
        }
    }
    else
    {
        tfvp = ctor(p, iF);

        // Check if constraint type override and store patchType if so
        if (patchTypeCstrIter != patchConstructorTable_().end())
        {
            tfvp.ref().patchType() = actualPatchType;
        }
    }

    if (patchFieldType == calculatedFvPatchField<Type>::typeName)
    {
        tfvp->setCalculatedMode();
    }

    return tfvp;
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
{
    const word patchFieldType(dict.lookup("type"));

    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTable_().find(patchFieldType);

    bool isGeneric = false;
    if (cstrIter == dictionaryConstructorTable_().end())
    {
        isGeneric = true;
        if (!disallowGenericFvPatchField)
        {
            cstrIter = dictionaryConstructorTable_().find("generic");
        }

        if (cstrIter == dictionaryConstructorTable_().end())
        {
            OSstream& s(FatalIOErrorInFunction(dict));
            s   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl
                << "Valid options are:" << nl << nl;
            printCtorTableKeys(s, dictionaryConstructorTable_());
            s   << exit(FatalIOError);
        }
    }

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter
            = dictionaryConstructorTable_().find(p.type());

        if (patchTypeCstrIter != dictionaryConstructorTable_().end() &&
            patchTypeCstrIter->second != cstrIter->second)
        {
            tmp<fvPatchField<Type>> tpf(cstrIter->second(p, iF, dict));
            typename dictionaryConstructorTable::iterator constraintTypeCstrIter
                = dictionaryConstructorTable_().find(tpf->constraintType());

            if (patchTypeCstrIter == constraintTypeCstrIter)
            {
                return tpf;
            }

            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    tmp<fvPatchField<Type>> tpf(cstrIter->second(p, iF, dict));

    if (!isGeneric)
    {
        // Warn about deprecated names
        wordList deprecatedPatchFields
        ({
            // List any deprecated patch names here
            // Currently empty
        });

        // Warn about deprecated compressible namespace
        word compNS = "compressible::";
        if
        (
            deprecatedPatchFields.found(patchFieldType)
         || patchFieldType.substr(0, compNS.size()) == compNS
        )
        {
            DeprecationWarningInFunction
            (
                patchFieldType,
                "boundary condition",
                40300,
                "Please use '" + tpf().type() + "' instead."
            );
        }
    }

    return tpf;
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvPatchField<Type>" << endl;
    }

    const auto ctor =
        ctorTableLookup
        (
            "patchField type", patchMapperConstructorTable_(), ptf.type()
        );
    return ctor(ptf, p, iF, pfMapper);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::fvPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    tmp<fvPatchField<Type>> tfvpf;

    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTable_().find(p.type());

    if (patchTypeCstrIter != patchConstructorTable_().end())
    {
        tfvpf =
            patchTypeCstrIter->second
            (
                p,
                DimensionedField<Type, volMesh>::null()
            );
    }
    else
    {
        tfvpf =
            tmp<fvPatchField<Type>>
            (
                new calculatedFvPatchField<Type>
                (
                    p,
                    DimensionedField<Type, volMesh>::null()
                )
            );
    }
    tfvpf->setCalculatedMode();
    return tfvpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
