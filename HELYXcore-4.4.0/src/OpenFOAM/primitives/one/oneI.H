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

#include "primitives/one/one.H"
#include "primitives/Scalar/scalar/scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class arg2>
class innerProduct<one, arg2>
{
public:

    typedef arg2 type;
};

inline const one& operator*(const one& o, const one&)
{
    return o;
}

template<class Type>
inline const Type& operator*(const Type& t, const one&)
{
    return t;
}

template<class Type>
inline const Type& operator*(const one&, const Type& t)
{
    return t;
}

template<class Type>
inline const Type& operator&(const one&, const Type& t)
{
    return t;
}

inline const one& operator/(const one& o, const one&)
{
    return o;
}

template<class Type>
inline Type operator/(const one&, const Type& t)
{
    return scalar(1)/t;
}

template<class Type>
inline const Type& operator/(const Type& t, const one&)
{
    return t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
