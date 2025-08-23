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
    Vector of complex numbers.

\*---------------------------------------------------------------------------*/

#include "primitives/Vector/complexVector/complexVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::complexVector::vsType::typeName = "complexVector";

template<>
const char* const Foam::complexVector::vsType::componentNames[] =
{
    "x", "y", "z"
};

template<>
const Foam::complexVector Foam::complexVector::vsType::zero
(
    complexVector::uniform(complex(0, 0))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::one
(
    complexVector::uniform(complex(1, 1))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::max
(
    complexVector::uniform(complex(VGREAT, VGREAT))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::min
(
    complexVector::uniform(complex(-VGREAT, -VGREAT))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMax
(
    complexVector::uniform(complex(ROOTVGREAT, ROOTVGREAT))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMin
(
    complexVector::uniform(complex(-ROOTVGREAT, -ROOTVGREAT))
);


// ************************************************************************* //
