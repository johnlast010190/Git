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

\*---------------------------------------------------------------------------*/

#include "primitives/zero/zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline const Type& operator+(const Type& t, const zero&)
{
    return t;
}

template<class Type>
inline const Type& operator+(const zero&, const Type& t)
{
    return t;
}

template<class Type>
inline const Type& operator-(const Type& t, const zero&)
{
    return t;
}

template<class Type>
inline Type operator-(const zero&, const Type& t)
{
    return -t;
}

template<class Type>
inline zero operator*(const Type& t, const zero&)
{
    return Zero;
}

template<class Type>
inline zero operator*(const zero&, const Type& t)
{
    return Zero;
}

template<class Type>
inline zero operator/(const zero&, const Type& t)
{
    return Zero;
}

template<typename Type>
Type operator/(const Type& t, const zero&)
{
    static_assert(!std::is_same_v<Type, Type>, "Division by zero!");
    return t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
