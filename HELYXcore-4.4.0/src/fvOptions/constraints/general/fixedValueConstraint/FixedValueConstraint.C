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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "FixedValueConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::FixedValueConstraint<Type>::FixedValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::FixedValueConstraint<Type>::sourceFields
(
    wordList& fieldNames
)
{
    const dictionary& fieldValuesDict = coeffs_.subDict("fieldValues");

    fieldNames.setSize(fieldValuesDict.size());
    fieldValues_.setSize(fieldNames.size());

    label i = 0;
    forAllConstIter(dictionary, fieldValuesDict, iter)
    {
        fieldNames[i] = iter().keyword();
        fieldValuesDict.lookup(iter().keyword()) >> fieldValues_[i];
        i++;
    }
}


template<class Type>
void Foam::fv::FixedValueConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInformation
        << "FixedValueConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    eqn.setValues(cells_, List<Type>(cells_.size(), fieldValues_[fieldi]));
}

template<class Type>
void Foam::fv::FixedValueConstraint<Type>::constrain
(
    fvBlockMatrix<Type>& eqn,
    const label fieldi
)
{
    NotImplemented;
}

// ************************************************************************* //
