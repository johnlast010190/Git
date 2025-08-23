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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::PvdFieldsType

Description
    Enum with the possible ways that the fields for the PVD export can be chosen

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class PvdFieldsType
    {
    public:
        enum Value
        {
            ALL,
            FROM_LIST,
            FROM_SCENE
        };

        PvdFieldsType(){ value = FROM_SCENE; }
        explicit PvdFieldsType(Value aValue) {value = aValue;}
        explicit PvdFieldsType(const word& aType) : value(fromWord(aType)) {}
        PvdFieldsType(const PvdFieldsType& type) {value = type.value;}

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
        bool operator==(const PvdFieldsType& i) const { return value == i.value; }
        bool operator==(Value i) const { return value == i; }
        bool operator!=(const PvdFieldsType& i) const { return value != i.value; }

        PvdFieldsType& operator=(const PvdFieldsType& i) = default;
        PvdFieldsType& operator=(const Value& i) {value = i; return *this;}

        friend Ostream& operator<<(Ostream& os, const PvdFieldsType& PvdFieldsType)
        {
            os << PvdFieldsType.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, PvdFieldsType& m)
        {
            word w;
            is >> w;
            m.value = fromWord(w);
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            const static HashTable<Value, word> dictPvdFieldsTypes_
                {
                    {"all", ALL},
                    {"fromlist", FROM_LIST},
                    {"fromscene", FROM_SCENE}
                };
            return dictPvdFieldsTypes_;
        }
        Value value;

        static Value fromWord(word w)
        {
            std::transform(w.begin(), w.end(), w.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (getTable().found(w))
            {
                return getTable()[w];
            }
            else
            {
                Warning
                    << "Unsupported PVD field selection type:"
                    << nl << indent << w << nl
                    << R"(Supported formats are "all", "fromList" and "fromScene")"
                    << endl;
                return ALL;
            }
        }
    } ;
} // End namespace Foam
