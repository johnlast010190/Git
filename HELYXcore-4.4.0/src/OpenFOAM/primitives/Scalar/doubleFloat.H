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
    (c) 2011 Symscape
    (c) 2014 blueCAPE Lda
    (c) 2016 blueCAPE Lda
    (c) 2025 Engys Ltd

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.



\*---------------------------------------------------------------------------*/

#ifndef doubleFloat_H
#define doubleFloat_H

#include "primitives/ints/label/label.H"
#include "primitives/VectorSpace/products.H"

#if (defined( WIN32 ) || defined( WIN64 )) && defined(__STRICT_ANSI__)
#undef __STRICT_ANSI__
#define REDEFINE_SA__
#endif

#include <cmath>

#if defined( WIN32 ) || defined( WIN64 )
#ifdef REDEFINE_SA__
#define __STRICT_ANSI__ 1
#endif

#define j0f _j0
#define j1f _j1
#define jnf _jn
#define y0f _y0
#define y1f _y1
#define ynf _yn
#endif

#ifdef HELYX_ALTERNATIVE_LIBM
#include "global/libmFunctions/libmFunctions.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline bool equal(const T& s1, const T& s2)
{
    return s1 == s2;
}


#ifdef HELYX_ALTERNATIVE_LIBM
#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return libmFunctions::pow(double(s), double(e)); \
}
#else
#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::pow(double(s), double(e));           \
}
#endif


MAXMINPOW(double, double, double)
MAXMINPOW(double, double, float)
MAXMINPOW(double, float, double)
MAXMINPOW(float, float, float)
MAXMINPOW(double, double, int)
MAXMINPOW(double, int, double)
MAXMINPOW(double, double, long)
MAXMINPOW(double, long, double)
MAXMINPOW(float, float, int)
MAXMINPOW(float, int, float)
MAXMINPOW(float, float, long)
MAXMINPOW(float, long, float)

#undef MAXMINPOW


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
