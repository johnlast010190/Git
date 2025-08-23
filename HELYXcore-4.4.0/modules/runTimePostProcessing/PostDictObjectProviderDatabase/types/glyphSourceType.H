/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::GlyphSourceType

Description
    Enum with the possible source types for the vectors

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"
#include "engysGlyphs.h"

namespace Foam::functionObjects::runTimeVis
{
    class GlyphSourceType
    {
    public:
        enum Value
        {
            ALL_POINTS = engysGlyphs::GetDistributionCodeForAllPoints(),
            UNIFORM_BOUNDS_BASED = engysGlyphs::GetDistributionCodeForUniform(),
            EVERY_NTH_POINT = engysGlyphs::GetDistributionCodeForEveryNthPoint()
        };

        GlyphSourceType(){ value = ALL_POINTS; }
        explicit GlyphSourceType(const word& aSource) : value(getTable()[aSource]) {}

        [[nodiscard]] int getValue() const {return value;}

        [[nodiscard]] word getName() const
        {
            for (auto coordType = getTable().begin(); coordType != getTable().end(); coordType++)
            {
                if (coordType.object() == value)
                {
                    return coordType.key();
                }
            }
            return "unknown";
        }

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const GlyphSourceType& f) const { return value == f.value; }
        bool operator!=(const GlyphSourceType& f) const { return value != f.value; }

        friend Ostream& operator<<(Ostream& os, const GlyphSourceType& type)
        {
            os << type.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, GlyphSourceType& m)
        {
            word w;
            is >> w;
            m.value = getTable()[w];
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> sources_
                {
                    {"allPoints", ALL_POINTS},
                    {"uniformBoundsBased", UNIFORM_BOUNDS_BASED},
                    {"everyNthPoint", EVERY_NTH_POINT},
                };
            return sources_;
        }
        Value value;
    } ;
} // End namespace
