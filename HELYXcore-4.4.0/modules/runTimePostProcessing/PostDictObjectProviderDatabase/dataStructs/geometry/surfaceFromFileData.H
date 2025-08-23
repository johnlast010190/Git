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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::SurfaceFromFileData

Description
    Contains the dict data necessary for the surfaces read from files, like STLs

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "postDict/postDictKeys.H"
#include "Utils/Utils.H"
#include "hash/hasher.H"
#include "db/dictionary/dictionary.H"
#include "primitives/strings/fileName/fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations

namespace Foam::functionObjects::runTimeVis
{

struct SurfaceFromFileData {
    string caseFolder;
    fileName fileNameWithExtension;
    string solidSubstring;
    bool mergeCoplanar;

    void readDict(const dictionary& itemDict, const dictionary& geometryDict, string caseFolder_)
    {
        fileNameWithExtension = Utils::getSubdictName(geometryDict);
        caseFolder = std::move(caseFolder_);

        string solidName = itemDict.lookupOrDefault(surfaceKeys::NAME_KEY, Utils::getSubdictName(itemDict));
        size_t surfaceFileNameSize = fileNameWithExtension.lessExt().length() + 1;
        if (solidName.length() > surfaceFileNameSize)
        {
            solidSubstring = solidName.substr(surfaceFileNameSize);
        }
        else
        {
            solidSubstring = "";
        }
        mergeCoplanar = itemDict.lookupOrDefault<bool>(surfaceKeys::MERGE_COPLANAR_KEY, false);
    }

    [[nodiscard]] bool readAllSurfaces() const
    {
        return solidSubstring.empty();
    }

    bool operator==(const SurfaceFromFileData& other) const
    {
        return caseFolder == other.caseFolder &&
               fileNameWithExtension == other.fileNameWithExtension &&
               solidSubstring == other.solidSubstring &&
               mergeCoplanar == other.mergeCoplanar;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, caseFolder);
        hasher::hash_combine(hash, fileNameWithExtension);
        hasher::hash_combine(hash, solidSubstring);
        hasher::hash_combine(hash, mergeCoplanar);
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
