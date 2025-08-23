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
    (c) 2020 Engys Ltd.

Typedef
    Foam::symmTransformField

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    symmTransformField.C

\*---------------------------------------------------------------------------*/

#ifndef symmTransformField_H
#define symmTransformField_H

#include "primitives/transform/symmTransform.H"
#include "fields/Fields/symmTensorField/symmTensorField.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void transform(Field<Type>&, const symmTensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const symmTensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const symmTensorField&, const tmp<Field<Type>>&);

template<class Type>
tmp<Field<Type>> transform(const tmp<symmTensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform
(
    const tmp<symmTensorField>&,
    const tmp<Field<Type>>&
);


template<class Type>
void transform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensor&, const tmp<Field<Type>>&);


template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const symmTensorField&
);

template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tmp<symmTensorField>&
);


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const symmTensorField&
);

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>&
);


template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const symmTensorField&
);

template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const tmp<symmTensorField>&
);


//- Return the mean of the original value and its value reflected in the
//  plane with normal n (need not be unit normal)
template<class Type>
tmp<Field<Type>>
symmetryValue(const Field<Type>&, const vectorField& n);

template<>
tmp<vectorField>
symmetryValue<vector>(const vectorField&, const vectorField& n);

template<>
tmp<scalarField>
symmetryValue<scalar>(const scalarField&, const vectorField& n);

template<class Type>
tmp<Field<Type>>
symmetryValue(const tmp<Field<Type>>&, const vectorField& n);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/Fields/symmTransformField/symmTransformField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
