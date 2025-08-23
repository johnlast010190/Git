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
    (c) 2017-2021 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

Typedef
    Foam::scalar

Description
    Single floating point number identical to float or double depending on
    whether HELYX_SP or HELYX_DP is defined.

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "primitives/Scalar/floatScalar/floatScalar.H"
#include "primitives/Scalar/doubleScalar/doubleScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(HELYX_SP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;

    static const scalar GREAT = floatScalarGREAT;
    static const scalar VGREAT = floatScalarVGREAT;
    static const scalar ROOTVGREAT = floatScalarROOTVGREAT;
    static const scalar SMALL = floatScalarSMALL;
    static const scalar ROOTSMALL = floatScalarROOTSMALL;
    static const scalar VSMALL = floatScalarVSMALL;
    static const scalar ROOTVSMALL = floatScalarROOTVSMALL;
    static const scalar NaN = floatScalarNaN;

    scalar readScalar(Istream& is);

    //- Read scalar from stream if present or return default value
    scalar readScalarOrDefault(Istream& is, const scalar defaultValue);
}

#elif defined(HELYX_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;

    static const scalar GREAT = doubleScalarGREAT;
    static const scalar VGREAT = doubleScalarVGREAT;
    static const scalar ROOTVGREAT = doubleScalarROOTVGREAT;
    static const scalar SMALL = doubleScalarSMALL;
    static const scalar ROOTSMALL = doubleScalarROOTSMALL;
    static const scalar VSMALL = doubleScalarVSMALL;
    static const scalar ROOTVSMALL = doubleScalarROOTVSMALL;
    static const scalar NaN = doubleScalarNaN;

    //- Read a single scalar from an input stream
    scalar readScalar(Istream& is);

    //- Read scalar from stream if present or return default value
    scalar readScalarOrDefault(Istream& is, const scalar defaultValue);

    //- Compute the power of the number x to the integer e
    inline scalar integerPow(const scalar x, const label e);

    //- Compute the power of the number x to the reciprocal integer 1/e
    inline scalar integerRoot(const scalar x, const label e);
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Additional transcendental functions

namespace Foam
{
    //- Inverse normalized incomplete gamma function
    scalar invIncGamma(const scalar a, const scalar P);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/Scalar/scalar/scalarI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
