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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceType

Description
    Enum with the possible surface types

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class SurfaceType
    {
    public:
        enum Value
        {
            FROM_FILE = 1,
            BOX,
            SPHERE,
            CYLINDER,
            PLANE,
            RING,
            FEATURE_LINE,
            UNKNOWN
        };

        SurfaceType() {value = FROM_FILE;};
        explicit SurfaceType(const word& aSurfaceType) { value = surfaceTypes_[aSurfaceType]; }
        explicit SurfaceType(const Value& value) { this->value = value; }

        [[nodiscard]] Value getValue() const {return value;};
        void setValue(Value v) {value = v;}

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const SurfaceType& l) const { return value == l.value; }
        bool operator!=(const SurfaceType& l) const { return value != l.value; }

        friend Ostream& operator<<(Ostream& os, const SurfaceType& surfaceType)
        {
            os << surfaceType.getValue();
            return os;
        };

        friend Istream& operator>>(Istream& is, SurfaceType& m)
        {
            word w;
            is.read(w);
            m.value = m.surfaceTypes_[w];
            return is;
        }

    private:
        HashTable<Value, word> surfaceTypes_
        {
            {"triSurfaceMesh", FROM_FILE},
            {"searchableBox", BOX},
            {"searchableSphere", SPHERE},
            {"searchableCylinder", CYLINDER},
            {"searchablePlane", PLANE},
            {"searchableRing", RING},
            {"featureline", FEATURE_LINE},
            {"unknown", UNKNOWN}
        };
        Value value;
    };
} // End namespace Foam
