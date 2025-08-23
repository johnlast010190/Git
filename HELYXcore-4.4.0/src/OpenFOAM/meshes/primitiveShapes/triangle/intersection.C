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

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveShapes/triangle/intersection.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::intersection::planarTol_ = 0.2;

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::intersection::direction,
        2
    >::names[] =
    {
        "vector",
        "contactSphere"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::intersection::algorithm,
        3
    >::names[] =
    {
        "fullRay",
        "halfRay",
        "visible"
    };
}

const Foam::NamedEnum<Foam::intersection::direction, 2>
Foam::intersection::directionNames_;

const Foam::NamedEnum<Foam::intersection::algorithm, 3>
Foam::intersection::algorithmNames_;


// ************************************************************************* //
