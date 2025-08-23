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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "ExplicitSetDataEntryValue.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::ExplicitSetDataEntryValue<Type>::setFieldData
(
    const dictionary& dict, wordList& fieldNames
)
{
    fieldNames.setSize(dict.toc().size());
    sourceData_.setSize(fieldNames.size());

    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        fieldNames[i] = iter().keyword();
        sourceData_.set
        (
            i,
            Function1<Type>::New
            (
                fieldNames[i],
                dict
            )
        );

        i++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::ExplicitSetDataEntryValue<Type>::ExplicitSetDataEntryValue
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    sourceData_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::ExplicitSetDataEntryValue<Type>::sourceFields(wordList& fieldNames)
{
    setFieldData(coeffs_.subDict("sourceData"), fieldNames);
}


template<class Type>
void Foam::fv::ExplicitSetDataEntryValue<Type>::setValue
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "ExplicitSetDataEntryValue<"<< pTraits<Type>::typeName
            << ">::setValue for source " << name_ << endl;
    }

    List<Type> values(cells_.size());

    UIndirectList<Type>(values, cells_)
        = sourceData_[fieldI].value(eqn.psi().mesh().time().value());

    eqn.setValues(cells_, values);
}


// ************************************************************************* //
