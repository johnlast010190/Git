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
    (c) 2025 Engys Ltd.

Description
    Optional binding of libm functions loaded from an alternative library than
    the one that was compiled against

\*---------------------------------------------------------------------------*/

#ifndef libmFunctions_H
#define libmFunctions_H

#ifdef HELYX_ALTERNATIVE_LIBM

#include <string>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Place all the static data in a struct so it is not repeated locally in each
// translation unit

struct libmFunctions
{
    // Counter for initialisation in each translation unit
    static int initCounter;

    static void* libHandle;

    template <class T>
    static T* libFunc(const std::string& fName, T* funcPtr);

    static double (*cbrt)(const double s);
    static double (*exp)(const double s);
    static double (*log)(const double s);
    static double (*log10)(const double s);
    static double (*sin)(const double s);
    static double (*cos)(const double s);
    static double (*tan)(const double s);
    static double (*asin)(const double s);
    static double (*acos)(const double s);
    static double (*atan)(const double s);
    static double (*sinh)(const double s);
    static double (*cosh)(const double s);
    static double (*tanh)(const double s);
    static double (*asinh)(const double s);
    static double (*acosh)(const double s);
    static double (*atanh)(const double s);
    static double (*erf)(const double s);
    static double (*erfc)(const double s);
    static double (*lgamma)(const double s);
    static double (*tgamma)(const double s);
    static double (*j0)(const double s);
    static double (*j1)(const double s);
    static double (*y0)(const double s);
    static double (*y1)(const double s);
    static double (*pow)(const double s, const double e);
    static double (*hypot)(const double s, const double e);
    static double (*atan2)(const double s, const double e);
    static double (*jn)(const int n, const double s);
    static double (*yn)(const int n, const double s);

    static float (*cbrtf)(const float s);
    static float (*expf)(const float s);
    static float (*logf)(const float s);
    static float (*log10f)(const float s);
    static float (*sinf)(const float s);
    static float (*cosf)(const float s);
    static float (*tanf)(const float s);
    static float (*asinf)(const float s);
    static float (*acosf)(const float s);
    static float (*atanf)(const float s);
    static float (*sinhf)(const float s);
    static float (*coshf)(const float s);
    static float (*tanhf)(const float s);
    static float (*asinhf)(const float s);
    static float (*acoshf)(const float s);
    static float (*atanhf)(const float s);
    static float (*erff)(const float s);
    static float (*erfcf)(const float s);
    static float (*lgammaf)(const float s);
    static float (*tgammaf)(const float s);
    static float (*j0f)(const float s);
    static float (*j1f)(const float s);
    static float (*y0f)(const float s);
    static float (*y1f)(const float s);
    static float (*powf)(const float s, const float e);
    static float (*hypotf)(const float s, const float e);
    static float (*atan2f)(const float s, const float e);
    static float (*jnf)(const int n, const float s);
    static float (*ynf)(const int n, const float s);

    static void init();
    static void cleanup();
};

// Initialiser object local to each translation unit. Based on 'Nifty counter'
// trick to ensure that libmFunctions is initialised before it is used in any TU
static struct libInitialiser
{
    libInitialiser();
    ~libInitialiser();
} li;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // HELYX_ALTERNATIVE_LIBM

#endif

// ************************************************************************* //
