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
    Foam::functionObjects::runTimeVis::SurfaceFormat

Description
    Enum with the possible surface formats

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/word/word.H"
#include "db/IOstreams/IOstreams.H"

namespace Foam::functionObjects::runTimeVis
{
    class SurfaceFormat
    {
    public:
        enum Value
        {
            STL_ASCII,
            STL_BINARY,
            STL_BINARY_ENGYS,
            STL_GZ_ASCII,
            OBJ,
            BTS,
            NAS,
            NAS_GZ,
            EMESH,
            STP,
            UNKNOWN
        };

        SurfaceFormat();
        explicit SurfaceFormat(const word& aSurfaceFormat);

        static SurfaceFormat fromNameWithExtension(const word& surfaceName);

        static bool hasValidExtension(const word &surfaceName);
        static bool hasMatchingNameAndValidExtension(const word& expectedName, const word& nameWithExtension);

        [[nodiscard]] Value getValue() const {return format.value;};
        [[nodiscard]] word getExtension() const {return format.extension;};

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const SurfaceFormat& l) const { return format.value == l.format.value; }
        bool operator!=(const SurfaceFormat& l) const { return format.value != l.format.value; }

        friend Ostream& operator<<(Ostream& os, const SurfaceFormat& SurfaceFormat)
        {
            os << SurfaceFormat.format.key;
            return os;
        };

        friend Istream& operator>>(Istream& is, SurfaceFormat& m)
        {
            word w;
            is.read(w);
            m.format = findFormatFromKey(w);
            return is;
        }

    private:
        struct SurfaceFormat_
        {
            Value value;
            const char* key;
            const char* extension;
            SurfaceFormat_(Value value, const char* key, const char* extension) :
            value(value), key(key), extension(extension)
            {}
        };

        SurfaceFormat_ format;

        static const std::vector<SurfaceFormat_>& getSurfaceFormats();
        static const SurfaceFormat_& getUnknownFormat();
        static const SurfaceFormat_& findFormatFromKey(const word& key);

        explicit SurfaceFormat(const SurfaceFormat_& format);
    };
} // End namespace Foam
