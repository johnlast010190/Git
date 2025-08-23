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

Description
    Vector2D of scalars.

\*---------------------------------------------------------------------------*/

#include "primitives/Vector2D/vector2DVector/vector2DVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::vector2DVector::vsType::typeName = "vector2DVector";

template<>
const char* const Foam::vector2DVector::vsType::componentNames[] = {"x", "y"};

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::vsType::zero
(
    vector2DVector::uniform(0)
);

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::one
(
    vector2DVector::uniform(1)
);

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::max
(
    vector2DVector::uniform(VGREAT)
);

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::min
(
    vector2DVector::uniform(-VGREAT)
);

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::rootMax
(
    vector2DVector::uniform(ROOTVGREAT)
);

template<>
const Foam::vector2DVector Foam::vector2DVector::vsType::rootMin
(
    vector2DVector::uniform(-ROOTVGREAT)
);


// ************************************************************************* //
