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
    Foam::functionObjects::runTimeVis::Location

Description
    Enum with the possible positions for the color legend

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class Location
    {
    public:
        enum Value
        {
            TOP,
            TOP_RIGHT,
            RIGHT,
            BOTTOM_RIGHT,
            BOTTOM,
            BOTTOM_LEFT,
            LEFT,
            TOP_LEFT,
            CUSTOM
        };

        Location(){ value = RIGHT; }

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const Location& l) const { return value == l.value; }
        bool operator!=(const Location& l) const { return value != l.value; }

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

        [[nodiscard]] bool IsLeft() const { return (value == TOP_LEFT) || (value == LEFT) || (value == BOTTOM_LEFT); }
        [[nodiscard]] bool IsCenterX() const { return (value == TOP) || (value == BOTTOM); }
        [[nodiscard]] bool IsRight() const { return (value == TOP_RIGHT) || (value == RIGHT) || (value == BOTTOM_RIGHT);}
        [[nodiscard]] bool IsBottom() const { return (value == BOTTOM_LEFT) || (value == BOTTOM) || (value == BOTTOM_RIGHT);}
        [[nodiscard]] bool IsCenterY() const { return (value == LEFT) || (value == RIGHT); }
        [[nodiscard]] bool IsTop() const { return (value == TOP_LEFT) || (value == TOP) || (value == TOP_RIGHT); }
        [[nodiscard]] bool IsCustom() const { return (value == CUSTOM); }

        friend Ostream& operator<<(Ostream& os, const Location& type)
        {
            os << type.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, Location& m)
        {
            word w;
            is >> w;
            m.value = getTable()[w];
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> locations_
                {
                    {"TOP",          TOP},
                    {"TOP_RIGHT",    TOP_RIGHT},
                    {"RIGHT",        RIGHT},
                    {"BOTTOM_RIGHT", BOTTOM_RIGHT},
                    {"BOTTOM",       BOTTOM},
                    {"BOTTOM_LEFT",  BOTTOM_LEFT},
                    {"LEFT",         LEFT},
                    {"TOP_LEFT",     TOP_LEFT},
                    {"custom",       CUSTOM}
                };
            return locations_;
        }
        Value value;
    } ;
} // End namespace
