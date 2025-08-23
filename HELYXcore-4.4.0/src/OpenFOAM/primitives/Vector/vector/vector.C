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
    Vector of scalars.

\*---------------------------------------------------------------------------*/

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::vector::vsType::typeName = "vector";

template<>
const char* const Foam::vector::vsType::componentNames[] = {"x", "y", "z"};

template<>
const Foam::vector Foam::vector::vsType::zero(vector::uniform(0));

template<>
const Foam::vector Foam::vector::vsType::one(vector::uniform(1));

template<>
const Foam::vector Foam::vector::vsType::max(vector::uniform(VGREAT));

template<>
const Foam::vector Foam::vector::vsType::min(vector::uniform(-VGREAT));

template<>
const Foam::vector Foam::vector::vsType::rootMax(vector::uniform(ROOTVGREAT));

template<>
const Foam::vector Foam::vector::vsType::rootMin(vector::uniform(-ROOTVGREAT));


// ************************************************************************* //
