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
    Foam::functionObjects::runTimeVis::CompressionLevelType

Description
    Enum with the possible export types for a scene

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class CompressionLevelType
    {
    public:
        enum Value
        {
            FASTEST,
            FAST,
            BALANCED,
            SMALL,
            SMALLEST,
            UNKNOWN
        };

        CompressionLevelType(){ value = BALANCED; }
        explicit CompressionLevelType(Value aValue) {value = aValue;}
        explicit CompressionLevelType(int aType) : value(fromInt(aType)) {}
        CompressionLevelType(const CompressionLevelType& type) {value = type.value;}

        [[nodiscard]] int getValue() const {return value;}

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const CompressionLevelType& i) const { return value == i.value; }
        bool operator==(Value i) const { return value == i; }
        bool operator!=(const CompressionLevelType& i) const { return value != i.value; }

        CompressionLevelType& operator=(const CompressionLevelType& i) = default;
        CompressionLevelType& operator=(const Value& i) {value = i; return *this;}

        friend Ostream& operator<<(Ostream& os, const CompressionLevelType& CompressionLevelType)
        {
            os << CompressionLevelType.getValue();
            return os;
        };

        friend Istream& operator>>(Istream& is, CompressionLevelType& m)
        {
            int w;
            is >> w;
            m.value = fromInt(w);
            return is;
        }

    private:
        Value value;

        static Value fromInt(int w)
        {
            if (w < UNKNOWN)
            {
                return Value(w);
            }
            else
            {
                Warning
                    << "The following unsupported compression level was specified:"
                    << nl << indent << w << nl
                    << R"(Supported levels are (fastest, biggest) 0, 1, 2, 3, 4 (slowest, smallest))"
                    << endl;
                return UNKNOWN;
            }
        }
    } ;
} // End namespace Foam
