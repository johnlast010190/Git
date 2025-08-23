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
    Foam::functionObjects::runTimeVis::FieldSamplingPvdFileData

Description
    Gets the pvd file path from the Function Object dictionary of a FieldSamplingFunctionObject

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/fileSource/pvdFileData.H"

#include "db/dictionary/dictionary.H"
#include "postDict/postDictKeys.H"
#include "Utils/Utils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct FieldSamplingPvdFileData : public PvdFileData
{
    FieldSamplingPvdFileData
        (
            const dictionary &functionObjectDict,
            const string &caseFolder
        ) : PvdFileData(readFileNameAndPath(functionObjectDict), caseFolder)
    {
    }

private:
    static string readFileNameAndPath(const dictionary &functionObjectDict)
    {
        auto fileName = functionObjectDict.lookup<string>(fieldSamplingFunctionObjectKeys::SAMPLING_FILE_KEY);
        string foName = Utils::getSubdictName(functionObjectDict);
        return "postProcessing" / foName / fileName;
    }
};

} // End namespace

// ************************************************************************* //
