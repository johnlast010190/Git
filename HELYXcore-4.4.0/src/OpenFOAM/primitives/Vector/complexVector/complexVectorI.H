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

Description
    complexVector specific part of 3D complexVector obtained from
    generic Vector.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline complexVector operator*(const complex& v1, const complexVector& v2)
{
    return complexVector
    (
        v1*v2.x(),
        v1*v2.y(),
        v1*v2.z()
    );
}


inline complexVector operator*(const complexVector& v2, const complex& v1)
{
    return complexVector
    (
        v1*v2.x(),
        v1*v2.y(),
        v1*v2.z()
    );
}


inline complexVector operator/(const complexVector& v1, const complex& v2)
{
    return complexVector
    (
        v1.x()/v2,
        v1.y()/v2,
        v1.z()/v2
    );
}


inline complexVector operator/(const complex& v1, const complexVector& v2)
{
    return complexVector
    (
        v1/v2.x(),
        v1/v2.y(),
        v1/v2.z()
    );
}


// complexVector dot product

inline complex operator&(const complexVector& v1, const complexVector& v2)
{
    return complex
    (
        v1.x()*v2.x().conjugate()
      + v1.y()*v2.y().conjugate()
      + v1.z()*v2.z().conjugate()
    );
}


// complexVector cross product

inline complexVector operator^(const complexVector& v1, const complexVector& v2)
{
    return complexVector
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// complexVector cross product

inline complexVector operator^(const vector& v1, const complexVector& v2)
{
    return complexVector
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
