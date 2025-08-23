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
    Foam::functionObjects::runTimeVis::GlyphType

Description
    Enum with the possible arrow types for the vectors

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"
#include "engysGlyphs.h"

namespace Foam::functionObjects::runTimeVis
{
    class GlyphType
    {
    public:
        enum Value
        {
            ARROW_3D = engysGlyphs::GetGlyphTypeCodeForArrow3D(),
            CONE = engysGlyphs::GetGlyphTypeCodeForCone(),
            ARROW_2D = engysGlyphs::GetGlyphTypeCodeForArrow2D(),
            FILLED_ARROW_2D = engysGlyphs::GetGlyphTypeCodeForFilledArrow2D(),
            TRIANGLE = engysGlyphs::GetGlyphTypeCodeForTriangle()
        };

        GlyphType(){ value = ARROW_3D; }

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
        bool operator==(const GlyphType& f) const { return value == f.value; }
        bool operator!=(const GlyphType& f) const { return value != f.value; }

        friend Ostream& operator<<(Ostream& os, const GlyphType& type)
        {
            os << type.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, GlyphType& m)
        {
            word w;
            is >> w;
            m.value = getTable()[w];
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> glyphs_
                {
                    {"arrow3d",       ARROW_3D},
                    {"cone",          CONE},
                    {"arrow2d",       ARROW_2D},
                    {"filledArrow2d", FILLED_ARROW_2D},
                    {"triangle",      TRIANGLE},
                };
            return glyphs_;
        }
        Value value;
    } ;
} // End namespace
