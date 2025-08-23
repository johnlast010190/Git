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

\*---------------------------------------------------------------------------*/

#include "primitives/Vector/labelVector/labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelVector::vsType::typeName = "labelVector";

template<>
const char* const Foam::labelVector::vsType::componentNames[] =
{
    "x", "y", "z"
};

template<>
const Foam::labelVector Foam::labelVector::vsType::zero
(
    labelVector::uniform(0)
);

template<>
const Foam::labelVector Foam::labelVector::vsType::one
(
    labelVector::uniform(1)
);

template<>
const Foam::labelVector Foam::labelVector::vsType::max
(
    labelVector::uniform(labelMax)
);

template<>
const Foam::labelVector Foam::labelVector::vsType::min
(
    labelVector::uniform(-labelMax)
);

template<>
const Foam::labelVector Foam::labelVector::vsType::rootMax
(
    labelVector::uniform(sqrt(scalar(labelMax)))
);

template<>
const Foam::labelVector Foam::labelVector::vsType::rootMin
(
    labelVector::uniform(-sqrt(scalar(labelMax)))
);


// ************************************************************************* //
