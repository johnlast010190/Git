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

#include "surfaceFormat.H"
#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
const std::vector<SurfaceFormat::SurfaceFormat_> &SurfaceFormat::getSurfaceFormats()
{
    static const std::vector<SurfaceFormat::SurfaceFormat_> surfaceFormats_
        {
            {STL_ASCII,        "ascii",  ".stl"},
            {STL_BINARY,       "binary", ".stl"},
            {STL_BINARY_ENGYS, "ebs",    ".ebs"},
            {STL_GZ_ASCII,     "gz",     ".stl.gz"},
            {OBJ,              "obj",    ".obj"},
            {BTS,              "bts",    ".bts"},
            {NAS,              "nas",    ".nas"},
            {NAS_GZ,           "nas",    ".nas.gz"},
            {EMESH,            "eMesh",  ".eMesh"},
            {STP,              "stp",    ".stp"},
            {STP,              "stp",    ".step"},
            {STP,              "step",   ".step"},
            {STP,              "step",   ".stp"}
        };
    return surfaceFormats_;
}

const SurfaceFormat::SurfaceFormat_& SurfaceFormat::getUnknownFormat()
{
    static const SurfaceFormat_ unknown(UNKNOWN, "", "");
    return unknown;
}

SurfaceFormat::SurfaceFormat() : format(getUnknownFormat())
{}

SurfaceFormat::SurfaceFormat(const word &aSurfaceFormat)
    : SurfaceFormat(findFormatFromKey(aSurfaceFormat))
{}

const SurfaceFormat::SurfaceFormat_ &SurfaceFormat::findFormatFromKey(const word &key)
{
    const auto &i = std::find_if(
        getSurfaceFormats().begin(), getSurfaceFormats().end(),
        [key](const SurfaceFormat_ &f) -> bool { return key == f.key; }
    );
    if (i == getSurfaceFormats().end())
    {
        FatalErrorInFunction << "Unknown surface file format: " << key
                             << abort(FatalError);
        return getUnknownFormat();
    }
    else
    {
        return *i;
    }
}

SurfaceFormat SurfaceFormat::fromNameWithExtension(const word &surfaceName)
{
    const auto &i = std::find_if(
        getSurfaceFormats().begin(), getSurfaceFormats().end(),
        [surfaceName](const SurfaceFormat_ &f) -> bool {
            return surfaceName.endsWith(f.extension);
        }
    );

    if (i == getSurfaceFormats().end())
    {
        FatalErrorInFunction << "Unknown surface file format for " << surfaceName
                             << abort(FatalError);
        return SurfaceFormat{getUnknownFormat()};
    }
    else
    {
        return SurfaceFormat{*i};
    }
}

bool SurfaceFormat::hasMatchingNameAndValidExtension(const word &expectedName, const word &nameWithExtension)
{
    if (!nameWithExtension.startsWith(expectedName))
    {
        return false;
    }
    word extension = nameWithExtension.substr(expectedName.size());

    const auto &i = std::find_if(
        getSurfaceFormats().begin(), getSurfaceFormats().end(),
        [extension](const SurfaceFormat_ &f) -> bool {
            return extension == f.extension;
        }
    );
    return i != getSurfaceFormats().end();
}

bool SurfaceFormat::hasValidExtension(const word &surfaceName)
{
    const auto &i = std::find_if(
        getSurfaceFormats().begin(), getSurfaceFormats().end(),
        [surfaceName](const SurfaceFormat_ &f) -> bool {
            return surfaceName.endsWith(f.extension);
        }
    );
    return i != getSurfaceFormats().end();
}

SurfaceFormat::SurfaceFormat(const SurfaceFormat_ &format)
    : format(format)
{}
} // End namespace Foam
