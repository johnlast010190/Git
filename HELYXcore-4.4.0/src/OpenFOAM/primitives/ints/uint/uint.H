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

Primitive
    uint

Description
    System uinteger

SourceFiles
    uintIO.C

\*---------------------------------------------------------------------------*/

#ifndef uint_H
#define uint_H

#include "primitives/ints/uint32/uint32.H"
#include "primitives/ints/uint64/uint64.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define MAXMIN(retType, type1, type2)              \
                                                   \
constexpr inline retType max(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 > s2)? s1: s2;                      \
}                                                  \
                                                   \
constexpr inline retType min(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 < s2)? s1: s2;                      \
}


MAXMIN(uint8_t, uint8_t, uint8_t)
MAXMIN(uint16_t, uint16_t, uint16_t)

MAXMIN(uint32_t, uint32_t, uint32_t)
MAXMIN(uint64_t, uint64_t, uint32_t)
MAXMIN(uint64_t, uint32_t, uint64_t)
MAXMIN(uint64_t, uint64_t, uint64_t)


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

unsigned int readUint(Istream&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
