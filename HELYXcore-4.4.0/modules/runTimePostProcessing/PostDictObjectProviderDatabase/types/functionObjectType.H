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
    Foam::functionObjects::runTimeVis::FunctionObjectType

Description
    Enum with the possible function object types that generate a file source, like a pvd

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class FunctionObjectType
    {
    public:
        enum Value
        {
            FIELD_SAMPLING,
            UNKNOWN
        };

        FunctionObjectType() {value = UNKNOWN;};

        [[nodiscard]] Value getValue() const {return value;};

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
        bool operator==(const FunctionObjectType& l) const { return value == l.value; }
        bool operator!=(const FunctionObjectType& l) const { return value != l.value; }

        friend Ostream& operator<<(Ostream& os, const FunctionObjectType& type)
        {
            os << type.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, FunctionObjectType& m)
        {
            word w;
            is >> w;
            m.value = getTable()[w];
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> functionTypes_
                {
                    {"fieldSampling", FIELD_SAMPLING},
                    {"unknown",       UNKNOWN}
                };
            return functionTypes_;
        }
        Value value;
    };
} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam
