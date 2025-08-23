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
    Foam::functionObjects::runTimeVis::ClipObjectData

Description
    Contains the dict data for a clip, but without the clipType subdict,
    which is held in one of the cuttingData instead.

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "hash/hasher.H"
#include "postDict/postDictKeys.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct ClipObjectData
{
    bool crinkle;
    bool insideOut;
    word referenceFrame;

    void readDict(const dictionary& objectDict)
    {
        crinkle = objectDict.lookup<bool>(clipObjectKeys::CRINKLE_KEY);
        insideOut = objectDict.lookupOrDefault<bool>(clipObjectKeys::INSIDE_OUT_KEY, false);
        referenceFrame = objectDict.lookupOrDefault<word>(cameraKeys::REFERENCE_FRAME_KEY, "");
    }

    bool operator==(const ClipObjectData& other) const
    {
        return crinkle == other.crinkle &&
             insideOut == other.insideOut &&
             referenceFrame == other.referenceFrame;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, crinkle);
        hasher::hash_combine(hash, insideOut);
        hasher::hash_combine(hash, referenceFrame);
        return hash;
    }
};

} // End namespace

// ************************************************************************* //
