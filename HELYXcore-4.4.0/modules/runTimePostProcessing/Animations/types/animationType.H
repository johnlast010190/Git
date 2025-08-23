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
    Foam::functionObjects::runTimeVis::AnimationType

Description
    Enum with the possible animation types

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::animations
{
    class AnimationType
    {
    public:
        enum Value
        {
            AVI,
            MP4,
            GIF,
            WEBP
        };

        AnimationType() {value = MP4;};
        explicit AnimationType(const word& aAnimationType) { value = animationTypes_[aAnimationType]; }

        [[nodiscard]] Value getValue() const {return value;};
        void setValue(Value v) {value = v;}

        [[nodiscard]] word getExtension() const {
            for (auto coordType = animationTypes_.begin(); coordType != animationTypes_.end(); coordType++)
            {
                if (coordType.object() == value)
                {
                    return coordType.key();
                }
            }
            return "unknown";
        };

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const AnimationType& l) const { return value == l.value; }
        bool operator!=(const AnimationType& l) const { return value != l.value; }

        friend Ostream& operator<<(Ostream& os, const AnimationType& surfaceType)
        {
            os << surfaceType.getValue();
            return os;
        };

        friend Istream& operator>>(Istream& is, AnimationType& m)
        {
            word w;
            is >> w;
            m.fromWord(w);
            return is;
        }

    private:
        HashTable<Value, word> animationTypes_
        {
            {"avi", AVI},
            {"mp4", MP4},
            {"gif", GIF},
            {"webp", WEBP}
        };
        Value value;

        void fromWord(word w)
        {
            std::transform(w.begin(), w.end(), w.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            value = animationTypes_[w];
        }
    };
} // End namespace Foam
