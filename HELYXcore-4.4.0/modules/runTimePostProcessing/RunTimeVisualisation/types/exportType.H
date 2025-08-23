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
    (c) 2022-2023 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ExportType

Description
    Enum with the possible export types for a scene

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class ExportType
    {
    public:
        enum Value
        {
            PNG,
            X3D,
            PVD,
            EDF,
            UNKNOWN
        };

        ExportType(){ value = PNG; }
        explicit ExportType(Value aValue) {value = aValue;}
        explicit ExportType(const word& aType) : value(fromWord(aType)) {}
        ExportType(const ExportType& type) {value = type.value;}

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
        bool operator==(const ExportType& i) const { return value == i.value; }
        bool operator==(Value i) const { return value == i; }
        bool operator!=(const ExportType& i) const { return value != i.value; }

        ExportType& operator=(const ExportType& i) = default;
        ExportType& operator=(const Value& i) {value = i; return *this;}

        friend Ostream& operator<<(Ostream& os, const ExportType& ExportType)
        {
            os << ExportType.getName();
            return os;
        };

        friend Istream& operator>>(Istream& is, ExportType& m)
        {
            word w;
            is >> w;
            m.value = fromWord(w);
            return is;
        }

    private:
        static const HashTable<Value, word>& getTable()
        {
            static HashTable<Value, word> dictExportTypes_
                {
                    {"png", PNG},
                    {"x3d", X3D},
                    {"pvd", PVD},
                    {"edf", EDF}
                };
            return dictExportTypes_;
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
                    << "The following unsupported export format was specified:"
                    << nl << indent << w << nl
                    << R"(Supported formats are "png", "x3d" and "pvd")"
                    << endl;
                return UNKNOWN;
            }
        }
    } ;
} // End namespace Foam
