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
    Foam::functionObjects::runTimeVis::BackfaceStyling

Description
    Enum with the possible backface styling for the items

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class BackfaceStyling
    {
    public:
        enum Value
        {
            FOLLOW_FRONTFACE,
            CULL_BACKFACE,
            CULL_FRONTFACE
        };

        BackfaceStyling(){ value = FOLLOW_FRONTFACE; }

        [[nodiscard]] Value getValue() const {return value;}

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
        bool operator==(const BackfaceStyling& b) const { return value == b.value; }
        bool operator!=(const BackfaceStyling& b) const { return value != b.value; }

        friend Ostream& operator<<(Ostream& os, const BackfaceStyling& type)
        {
            os << type.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, BackfaceStyling& m)
        {
            word w;
            is >> w;
            m.value = getTable()[w];
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> stylings_
                {
                    {"followFrontface", FOLLOW_FRONTFACE},
                    {"cullBackface",    CULL_BACKFACE},
                    {"cullFrontface",   CULL_FRONTFACE},
                };
            return stylings_;
        }
        Value value;
    } ;
} // End namespace
