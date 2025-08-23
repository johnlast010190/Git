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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FromFileLine3DSourceData

Description
    Contains the dict data necessary to define a line source read from a file for the line 3D

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "types/line3DSourceType.H"
#include "hash/hasher.H"
#include "baseClasses/id.H"
#include "itemDataSetProviders/itemDataSetProvider.H"

#include "db/dictionary/dictionary.H"

#include "engysLine3D.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct FromFileLine3DSourceData
{
    fileName filePath;

    explicit FromFileLine3DSourceData(const dictionary& line3DDict, const string& caseFolder)
    {
        readDict(line3DDict, caseFolder);
    }

    void readDict(const dictionary& line3DDict, const string& caseFolder)
    {
        auto name = line3DDict.lookup<string>(line3DObjectKeys::FILE_NAME_KEY);
        filePath = caseFolder / name;
    }

    bool operator==(const FromFileLine3DSourceData& other) const
    {
        return filePath == other.filePath;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, Line3DSourceType::FROM_FILE);
        hasher::hash_combine(hash, filePath);
    }
};

} // End namespace Foam

// ************************************************************************* //
