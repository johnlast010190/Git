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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::bjSlopeLimiter

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Class with limiter function which returns the slope limiter for the
    Barth-Jespersen differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

    Original source:
    T. Barth and D. Jespersen, 1989.
    "The design and application of upwind schemes on unstructured meshes."
    https://arc.aiaa.org/doi/abs/10.2514/6.1989-366

SourceFiles
    bjSlope.C

\*---------------------------------------------------------------------------*/

#ifndef bjSlope_H
#define bjSlope_H

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class bjSlopeLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class bjSlopeLimiter
:
    public LimiterFunc
{

public:

    bjSlopeLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        const scalar r =
            LimiterFunc::r(faceFlux, phiP, phiN, gradcP, gradcN, d);

        return max(min(min(1, 4/(r + 1 + VSMALL)), 4*r/(r + 1 + VSMALL)), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
