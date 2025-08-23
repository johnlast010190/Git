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
    (c) 2018 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

InNamespace
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformField.C

\*---------------------------------------------------------------------------*/

#ifndef transformField_H
#define transformField_H

#include "primitives/transform/transform.H"
#include "primitives/spatialVectorAlgebra/spatialTransform/spatialTransform.H"
#include "primitives/quaternion/quaternion.H"
#include "primitives/septernion/septernion.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "fields/Fields/tensorField/tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
List<Type> transform(const tensor&, const UList<Type>&);


template<class Type>
void transform(Field<Type>&, const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensorField&, const tmp<Field<Type>>&);

template<class Type>
tmp<Field<Type>> transform(const tmp<tensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tmp<tensorField>&, const tmp<Field<Type>>&);


template<class Type>
void transform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> transform(const tensor&, const tmp<Field<Type>>&);

template<class Type>
void invTransform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
void invTransform(Field<Type>&, const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> invTransform(const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> invTransform(const tensorField&, const tmp<Field<Type>>&);

template<class Type>
tmp<Field<Type>> invTransform(const tmp<tensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> invTransform(const tmp<tensorField>&, const tmp<Field<Type>>&);

template<class Type>
tmp<Field<Type>> invTransform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type>> invTransform(const tensor&, const tmp<Field<Type>>&);

template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const Field<Type2>&);

template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const tmp<Field<Type2>>&);


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tensorField&
);

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<tensorField>&
);


template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tensorField&
);

template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>&
);


//- Transform given vectorField of coordinates with the given spatialTransform
void transformPoints(vectorField&, const spatialTransform&, const vectorField&);

//- Transform given vectorField of coordinates with the given spatialTransform
tmp<vectorField> transformPoints(const spatialTransform&, const vectorField&);

//- Transform given tmp<vectorField> of coordinates
//  with the given spatialTransform
tmp<vectorField> transformPoints
(
    const spatialTransform&,
    const tmp<vectorField>&
);


//- Rotate given vectorField with the given quaternion
void transform(vectorField&, const quaternion&, const vectorField&);

//- Rotate given vectorField with the given quaternion
tmp<vectorField> transform(const quaternion&, const vectorField&);

//- Rotate given tmp<vectorField> with the given quaternion
tmp<vectorField> transform(const quaternion&, const tmp<vectorField>&);


//- Transform given vectorField of coordinates with the given septernion
void transformPoints(vectorField&, const septernion&, const vectorField&);

//- Transform given vectorField of coordinates with the given septernion
tmp<vectorField> transformPoints(const septernion&, const vectorField&);

//- Transform given tmp<vectorField> of coordinates with the given septernion
tmp<vectorField> transformPoints(const septernion&, const tmp<vectorField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/Fields/transformField/transformFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
