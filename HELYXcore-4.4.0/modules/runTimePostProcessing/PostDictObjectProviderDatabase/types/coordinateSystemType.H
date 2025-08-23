/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::CoordinateSystemType

Description
    Enum with the possible coordinate types

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class CoordinateSystemType
    {
    public:
        enum Value
        {
            CARTESIAN,
            CYLINDRICAL
        };

        CoordinateSystemType(){ value = CARTESIAN; }
        explicit CoordinateSystemType(Value aValue) {value = aValue;}
        explicit CoordinateSystemType(const word& aType) { value = dictCoordTypes_[aType]; }
        CoordinateSystemType(const CoordinateSystemType& i) { value = i.value; }

        [[nodiscard]] int getValue() const {return value;}

        [[nodiscard]] word getName() const
        {
            for (auto coordType = dictCoordTypes_.begin(); coordType != dictCoordTypes_.end(); coordType++)
            {
                if (coordType.object() == value)
                {
                    return coordType.key();
                }
            }
            return "unknown";
        }

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const CoordinateSystemType& i) const { return value == i.value; }
        bool operator!=(const CoordinateSystemType& i) const { return value != i.value; }

        CoordinateSystemType& operator=(const CoordinateSystemType& i) {value = i.value; return *this;}
        CoordinateSystemType& operator=(const Value& i) {value = i; return *this;}

        friend Ostream& operator<<(Ostream& os, const CoordinateSystemType& coordType)
        {
            os << coordType.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, CoordinateSystemType& m)
        {
            word w;
            is.read(w);
            m.value = m.dictCoordTypes_[w];
            return is;
        }

    private:
        HashTable<Value, word> dictCoordTypes_
        {
            {"cartesian"  , CARTESIAN  },
            {"cylindrical", CYLINDRICAL},
        };
        Value value;
    } ;
} // End namespace Foam
