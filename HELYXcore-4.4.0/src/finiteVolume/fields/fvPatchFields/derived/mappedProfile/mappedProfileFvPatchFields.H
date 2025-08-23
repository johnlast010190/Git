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
    (c) 2012 OpenFOAM Foundation
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef mappedProfileFvPatchFields_H
#define mappedProfileFvPatchFields_H

#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchField.H"
#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef mappedProfileFvPatchField<scalar> mappedProfileFvPatchScalarField;
typedef mappedProfileFvPatchField<vector> mappedProfileFvPatchVectorField;
//makePatchTypeFieldTypedefs(mappedProfile);

template<>
tmp<Field<scalar>> mappedProfileFvPatchField<scalar>::mapProfileToField
(
    const pointField&,
    word, word
);

template <>
tmp<Field<vector>> mappedProfileFvPatchField<vector>::mapProfileToField
(
    const pointField& targetPoints,
    word fieldName,
    word regionName
);

template<>
void Foam::mappedProfileFvPatchField<scalar>::setTypeValue
(
    DynamicList<scalar>& field,
    const List<label>& locField,
    const List<string>& splitted
);

template<>
void Foam::mappedProfileFvPatchField<scalar>::transformTypeValue
(
    scalar& data,
    const vector& global
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
