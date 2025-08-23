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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Engys Ltd

Typedef
    Foam::vector

Description
    A scalar version of the templated Vector

SourceFiles
    vector.C

\*---------------------------------------------------------------------------*/

#ifndef vector_H
#define vector_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Vector/Vector.H"
#include "primitives/contiguous/contiguous.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Vector<scalar> vector;
typedef Tuple2<vector, vector> vectorTuple;

//- Data associated with vector type are contiguous
template<>
constexpr bool contiguous<vector>() {return true;}


template<class Type>
class flux
:
    public innerProduct<vector, Type>
{};


template<>
class flux<scalar>
{
public:

    typedef scalar type;
};


// Declare specialisations of static members assigned in the .C file

template<>
const char* const Foam::vector::vsType::typeName;

template<>
const char* const Foam::vector::vsType::componentNames[];

template<>
const Foam::vector Foam::vector::vsType::zero;

template<>
const Foam::vector Foam::vector::vsType::one;

template<>
const Foam::vector Foam::vector::vsType::max;

template<>
const Foam::vector Foam::vector::vsType::min;

template<>
const Foam::vector Foam::vector::vsType::rootMax;

template<>
const Foam::vector Foam::vector::vsType::rootMin;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
