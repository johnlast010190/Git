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

\*---------------------------------------------------------------------------*/

#ifdef HELYX_ALTERNATIVE_LIBM

#include "libmFunctions.H"
#include "include/OSspecific.H"

#include <cmath>

namespace Foam
{

// * * * * * * * * * * * Static data initialisation  * * * * * * * * * * * * //

void* libmFunctions::libHandle = nullptr;
int libmFunctions::initCounter = 0;

double (*libmFunctions::cbrt)(const double s) = nullptr;
double (*libmFunctions::exp)(const double s) = nullptr;
double (*libmFunctions::log)(const double s) = nullptr;
double (*libmFunctions::log10)(const double s) = nullptr;
double (*libmFunctions::sin)(const double s) = nullptr;
double (*libmFunctions::cos)(const double s) = nullptr;
double (*libmFunctions::tan)(const double s) = nullptr;
double (*libmFunctions::asin)(const double s) = nullptr;
double (*libmFunctions::acos)(const double s) = nullptr;
double (*libmFunctions::atan)(const double s) = nullptr;
double (*libmFunctions::sinh)(const double s) = nullptr;
double (*libmFunctions::cosh)(const double s) = nullptr;
double (*libmFunctions::tanh)(const double s) = nullptr;
double (*libmFunctions::asinh)(const double s) = nullptr;
double (*libmFunctions::acosh)(const double s) = nullptr;
double (*libmFunctions::atanh)(const double s) = nullptr;
double (*libmFunctions::erf)(const double s) = nullptr;
double (*libmFunctions::erfc)(const double s) = nullptr;
double (*libmFunctions::lgamma)(const double s) = nullptr;
double (*libmFunctions::tgamma)(const double s) = nullptr;
double (*libmFunctions::j0)(const double s) = nullptr;
double (*libmFunctions::j1)(const double s) = nullptr;
double (*libmFunctions::y0)(const double s) = nullptr;
double (*libmFunctions::y1)(const double s) = nullptr;
double (*libmFunctions::pow)(const double s, const double e) = nullptr;
double (*libmFunctions::hypot)(const double s, const double e) = nullptr;
double (*libmFunctions::atan2)(const double s, const double e) = nullptr;
double (*libmFunctions::jn)(const int n, const double s) = nullptr;
double (*libmFunctions::yn)(const int n, const double s) = nullptr;

float (*libmFunctions::cbrtf)(const float s) = nullptr;
float (*libmFunctions::expf)(const float s) = nullptr;
float (*libmFunctions::logf)(const float s) = nullptr;
float (*libmFunctions::log10f)(const float s) = nullptr;
float (*libmFunctions::sinf)(const float s) = nullptr;
float (*libmFunctions::cosf)(const float s) = nullptr;
float (*libmFunctions::tanf)(const float s) = nullptr;
float (*libmFunctions::asinf)(const float s) = nullptr;
float (*libmFunctions::acosf)(const float s) = nullptr;
float (*libmFunctions::atanf)(const float s) = nullptr;
float (*libmFunctions::sinhf)(const float s) = nullptr;
float (*libmFunctions::coshf)(const float s) = nullptr;
float (*libmFunctions::tanhf)(const float s) = nullptr;
float (*libmFunctions::asinhf)(const float s) = nullptr;
float (*libmFunctions::acoshf)(const float s) = nullptr;
float (*libmFunctions::atanhf)(const float s) = nullptr;
float (*libmFunctions::erff)(const float s) = nullptr;
float (*libmFunctions::erfcf)(const float s) = nullptr;
float (*libmFunctions::lgammaf)(const float s) = nullptr;
float (*libmFunctions::tgammaf)(const float s) = nullptr;
float (*libmFunctions::j0f)(const float s) = nullptr;
float (*libmFunctions::j1f)(const float s) = nullptr;
float (*libmFunctions::y0f)(const float s) = nullptr;
float (*libmFunctions::y1f)(const float s) = nullptr;
float (*libmFunctions::powf)(const float s, const float e) = nullptr;
float (*libmFunctions::hypotf)(const float s, const float e) = nullptr;
float (*libmFunctions::atan2f)(const float s, const float e) = nullptr;
float (*libmFunctions::jnf)(const int n, const float s) = nullptr;
float (*libmFunctions::ynf)(const int n, const float s) = nullptr;


// * * * * * * * * * * * * Static class members  * * * * * * * * * * * * * * //

template <class T>
T* libmFunctions::libFunc(const std::string& fName, T* funcPtr)
{
    return reinterpret_cast<T*>(dlSym(libHandle, fName));
}

void libmFunctions::init()
{
    libHandle = dlOpen(HELYX_ALTERNATIVE_LIBM);

    cbrt = libmFunctions::libFunc("cbrt", ::cbrt);
    exp = libmFunctions::libFunc("exp", ::exp);
    log = libmFunctions::libFunc("log", ::log);
    log10 = libmFunctions::libFunc("log10", ::log10);
    sin = libmFunctions::libFunc("sin", ::sin);
    cos = libmFunctions::libFunc("cos", ::cos);
    tan = libmFunctions::libFunc("tan", ::tan);
    asin = libmFunctions::libFunc("asin", ::asin);
    acos = libmFunctions::libFunc("acos", ::acos);
    atan = libmFunctions::libFunc("atan", ::atan);
    sinh = libmFunctions::libFunc("sinh", ::sinh);
    cosh = libmFunctions::libFunc("cosh", ::cosh);
    tanh = libmFunctions::libFunc("tanh", ::tanh);
    asinh = libmFunctions::libFunc("asinh", ::asinh);
    acosh = libmFunctions::libFunc("acosh", ::acosh);
    atanh = libmFunctions::libFunc("atanh", ::atanh);
    erf = libmFunctions::libFunc("erf", ::erf);
    erfc = libmFunctions::libFunc("erfc", ::erfc);
    lgamma = libmFunctions::libFunc("lgamma", ::lgamma);
    tgamma = libmFunctions::libFunc("tgamma", ::tgamma);
    j0 = libmFunctions::libFunc("j0", ::j0);
    j1 = libmFunctions::libFunc("j1", ::j1);
    y0 = libmFunctions::libFunc("y0", ::y0);
    y1 = libmFunctions::libFunc("y1", ::y1);
    pow = libmFunctions::libFunc("pow", ::pow);
    hypot = libmFunctions::libFunc("hypot", ::hypot);
    atan2 = libmFunctions::libFunc("atan2", ::atan2);
    jn = libmFunctions::libFunc("jn", ::jn);
    yn = libmFunctions::libFunc("yn", ::yn);

    cbrtf = libmFunctions::libFunc("cbrtf", ::cbrtf);
    expf = libmFunctions::libFunc("expf", ::expf);
    logf = libmFunctions::libFunc("logf", ::logf);
    log10f = libmFunctions::libFunc("log10f", ::log10f);
    sinf = libmFunctions::libFunc("sinf", ::sinf);
    cosf = libmFunctions::libFunc("cosf", ::cosf);
    tanf = libmFunctions::libFunc("tanf", ::tanf);
    asinf = libmFunctions::libFunc("asinf", ::asinf);
    acosf = libmFunctions::libFunc("acosf", ::acosf);
    atanf = libmFunctions::libFunc("atanf", ::atanf);
    sinhf = libmFunctions::libFunc("sinhf", ::sinhf);
    coshf = libmFunctions::libFunc("coshf", ::coshf);
    tanhf = libmFunctions::libFunc("tanhf", ::tanhf);
    asinhf = libmFunctions::libFunc("asinhf", ::asinhf);
    acoshf = libmFunctions::libFunc("acoshf", ::acoshf);
    atanhf = libmFunctions::libFunc("atanhf", ::atanhf);
    erff = libmFunctions::libFunc("erff", ::erff);
    erfcf = libmFunctions::libFunc("erfcf", ::erfcf);
    lgammaf = libmFunctions::libFunc("lgammaf", ::lgammaf);
    tgammaf = libmFunctions::libFunc("tgammaf", ::tgammaf);
    j0f = libmFunctions::libFunc("j0f", ::j0f);
    j1f = libmFunctions::libFunc("j1f", ::j1f);
    y0f = libmFunctions::libFunc("y0f", ::y0f);
    y1f = libmFunctions::libFunc("y1f", ::y1f);
    powf = libmFunctions::libFunc("powf", ::powf);
    hypotf = libmFunctions::libFunc("hypotf", ::hypotf);
    atan2f = libmFunctions::libFunc("atan2f", ::atan2f);
    jnf = libmFunctions::libFunc("jnf", ::jnf);
    ynf = libmFunctions::libFunc("ynf", ::ynf);
}


void libmFunctions::cleanup()
{
    dlClose(libHandle);
}


// * * * * * * * * * * * * * * Constructor/Destructor  * * * * * * * * * * * //

libInitialiser::libInitialiser()
{
    // Initially, initCounter was zero-initialised
    if (!libmFunctions::initCounter)
    {
        libmFunctions::init();
    }
    libmFunctions::initCounter++;
}

libInitialiser::~libInitialiser()
{
    libmFunctions::initCounter--;
    if (!libmFunctions::initCounter)
    {
        libmFunctions::cleanup();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // HELYX_ALTERNATIVE_LIBM

// ************************************************************************* //
