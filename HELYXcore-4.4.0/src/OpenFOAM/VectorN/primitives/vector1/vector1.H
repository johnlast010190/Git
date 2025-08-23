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
    (c) 2010 Ivor Clifford

Type
    Vector1

Description
    Vector1 obtained from generic VectorN

SourceFiles
    Vector1.C

\*---------------------------------------------------------------------------*/

#ifndef Vector1_H
#define Vector1_H

#include "primitives/Scalar/scalar/scalar.H"
#include "VectorN/primitives/VectorN.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef VectorN<scalar, 1> vector1;

//- Specify data associated with vector1 type is contiguous
template<>
constexpr bool contiguous<vector1>() {return true;}

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- overide  VectorSpace functions that assume more than 1 component
inline scalar magSqr(const vector1& v1)
{
    scalar s = magSqr(v1.v_[0]);
    return s;
}

inline scalar cmptMax(const vector1& v1)
{
    return v1.v_[0];
}

inline scalar cmptMin(const vector1& v1)
{
    return v1.v_[0];
}

inline scalar cmptSum(const vector1& v1)
{
    return v1.v_[0];
}

inline scalar cmptProduct(const vector1& v1)
{
    return v1.v_[0];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// Include inline implementations
#include "VectorN/primitives/vector1/vector1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
