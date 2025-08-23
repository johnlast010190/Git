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

Class
    Foam::intersection

Description
    Foam::intersection

SourceFiles
   intersection.C

\*---------------------------------------------------------------------------*/

#ifndef intersection_H
#define intersection_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class intersection Declaration
\*---------------------------------------------------------------------------*/

class intersection
{
private:

    // Static data

        //- Relative tolerance for point in a plane.
        //  The tolerance is relative to the object size.
        //  Used to determine if a half-ray from a point in plane of triangle
        // intersects the triangle or not.
        static scalar planarTol_;

public:

    enum direction
    {
        VECTOR,
        CONTACT_SPHERE
    };

    enum algorithm
    {
        FULL_RAY,        // Intersecting with full ray
        HALF_RAY,        // Intersecting with half ray
        VISIBLE          // Intersecting with the visible side
    };

    // Static Member Functions

        //- Direction names
        static const NamedEnum<direction, 2> directionNames_;

        //- Projection algorithm names
        static const NamedEnum<algorithm, 3> algorithmNames_;

        //- Return planar tolerance
        static scalar planarTol()
        {
            return planarTol_;
        }

        //- Set the planar tolerance, returning the previous value
        static scalar setPlanarTol(const scalar t)
        {
            if (t < -VSMALL)
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }

            scalar oldTol = planarTol_;
            planarTol_ = t;

            return oldTol;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
