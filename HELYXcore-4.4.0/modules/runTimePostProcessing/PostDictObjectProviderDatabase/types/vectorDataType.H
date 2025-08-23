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
    Foam::functionObjects::runTimeVis::VectorDataType

Description
    Enum with the possible vector widget data sources

SourceFiles
    vectorDataType.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class VectorDataType
    {
    public:
        enum Value
        {
            BODY_ACCELERATION,
            GRAVITY,
            BODY_ACCELERATION_AND_GRAVITY,
            NUMBER_OF_VALUES
        };

        VectorDataType() { value = BODY_ACCELERATION;};
        explicit VectorDataType(const word& aVectorType) { value = valueFromWord(aVectorType); }

        [[nodiscard]] Value getValue() const {return value;};
        [[nodiscard]] string getTitle() const {return titleFromValue(value) + " (m/s²):";};

        [[nodiscard]] bool hasGravity() const {return value == GRAVITY || value == BODY_ACCELERATION_AND_GRAVITY; }
        [[nodiscard]] bool hasBodyAcceleration() const {return value == BODY_ACCELERATION || value == BODY_ACCELERATION_AND_GRAVITY; }

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(VectorDataType l) const { return value == l.value; }
        bool operator!=(VectorDataType l) const { return value != l.value; }

        friend Ostream& operator<<(Ostream& os, const VectorDataType& type)
        {
            os << titleFromValue(type.value);
            return os;
        };

        friend Istream& operator>>(Istream& is, VectorDataType& m)
        {
            word w;
            is >> w;
            m.value = valueFromWord(w);
            return is;
        }

    private:

        static Value valueFromWord(const word& aVectorType)
        {
            static const HashTable<VectorDataType::Value, word> vectorTypes =
                {
                    {"bodyAcceleration", VectorDataType::BODY_ACCELERATION},
                    {"gravity", VectorDataType::GRAVITY},
                    {"bodyAccelerationGravity", VectorDataType::BODY_ACCELERATION_AND_GRAVITY}
                };
            return vectorTypes[aVectorType];
        }

        static string titleFromValue(Value value)
        {
            static const string titles[VectorDataType::NUMBER_OF_VALUES] =
                {
                    "Body Acceleration",
                    "Gravity",
                    "Total Acceleration"
                };
            return titles[value];
        }

        Value value;
    };
} // End namespace
