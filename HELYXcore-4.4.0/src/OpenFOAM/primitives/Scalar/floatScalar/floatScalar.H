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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Typedef
    Foam::floatScalar

Description
    Float precision floating point scalar type.

SourceFiles
    floatScalar.C

\*---------------------------------------------------------------------------*/

#ifndef floatScalar_H
#define floatScalar_H

#include "primitives/Scalar/doubleFloat.H"
#include "primitives/direction/direction.H"
#include "primitives/strings/word/word.H"

#include <limits>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef float floatScalar;

// Largest and smallest scalar values allowed in certain parts of the code
// (6 is the number of significant figures in an
//  IEEE single precision number.  See limits.h or float.h)
constexpr floatScalar floatScalarGREAT = 1.0e+6;
constexpr floatScalar floatScalarVGREAT = 1.0e+37;
constexpr floatScalar floatScalarROOTVGREAT = 1.0e+18;
constexpr floatScalar floatScalarSMALL = 1.0e-6;
constexpr floatScalar floatScalarROOTSMALL = 1.0e-3;
constexpr floatScalar floatScalarVSMALL = 1.0e-37;
constexpr floatScalar floatScalarROOTVSMALL = 1.0e-18;

constexpr floatScalar floatScalarNaN =
    numeric_limits<floatScalar>::signaling_NaN();

//- Read whole of buf as a scalar. Return true if succesful.
inline bool readScalar(const char* buf, floatScalar& s)
{
    char* endPtr;
    s = strtof(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar floatScalar
#define ScalarVGREAT floatScalarVGREAT
#define ScalarVSMALL floatScalarVSMALL
#define ScalarROOTVGREAT floatScalarROOTVGREAT
#define ScalarROOTVSMALL floatScalarROOTVSMALL
#define readScalar readFloatScalar

inline Scalar mag(const Scalar s)
{
    return ::fabsf(s);
}

#define stdTransFunc(func)         \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func##f(s);           \
}

#ifdef HELYX_ALTERNATIVE_LIBM
#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return libmFunctions::func##f(s); \
}
#else
#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func##f(s);           \
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/Scalar/Scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Scalar hypot(const Scalar x, const Scalar y)
{
#ifdef HELYX_ALTERNATIVE_LIBM
    return libmFunctions::hypotf(x, y);
#else
    return ::hypotf(x, y);
#endif
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
#ifdef HELYX_ALTERNATIVE_LIBM
    return libmFunctions::atan2f(y, x);
#else
    return ::atan2f(y, x);
#endif
}

inline Scalar jn(const int n, const Scalar s)
{
#ifdef HELYX_ALTERNATIVE_LIBM
    return libmFunctions::jnf(n, s);
#else
    return ::jnf(n, s);
#endif
}

inline Scalar yn(const int n, const Scalar s)
{
#ifdef HELYX_ALTERNATIVE_LIBM
    return libmFunctions::ynf(n, s);
#else
    return ::ynf(n, s);
#endif
}

#undef Scalar
#undef ScalarVGREAT
#undef ScalarVSMALL
#undef ScalarROOTVGREAT
#undef ScalarROOTVSMALL
#undef readScalar
#undef stdTransFunc
#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
