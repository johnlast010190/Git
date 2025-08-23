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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/ints/uint32/uint32.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const uint32_t Foam::pTraits<uint32_t>::zero = 0;
const uint32_t Foam::pTraits<uint32_t>::one = 1;
const uint32_t Foam::pTraits<uint32_t>::min = 0;
const uint32_t Foam::pTraits<uint32_t>::max = UINT32_MAX;
const uint32_t Foam::pTraits<uint32_t>::rootMin = 0;
const uint32_t Foam::pTraits<uint32_t>::rootMax = pTraits<uint32_t>::max;

const char* const Foam::pTraits<uint32_t>::componentNames[] = { "" };

Foam::pTraits<uint32_t>::pTraits(const uint32_t& p)
:
    p_(p)
{}

Foam::pTraits<uint32_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
