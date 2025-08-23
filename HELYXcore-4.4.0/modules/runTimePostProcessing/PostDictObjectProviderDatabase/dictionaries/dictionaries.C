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

\*---------------------------------------------------------------------------*/

#include "dictionaries.H"

#include "Utils/ParallelUtils.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "postDict/postDictKeys.H"

static const char *const HELYXHESHMESHDICT_NAME = "helyxHexMeshDict";
static const char *const GEOMETRYDICT_NAME = "geometryDict";
static const char *const BLOCKMESHDICT_NAME = "blockMeshDict";
static const char *const CONTROL_DICT_NAME = "controlDict";
static const char *const RTPP_DICT_NAME = "runtimeVisualisationDict";
static const char *const MESH_OBJECTS_DICT_NAME = "meshObjects";

namespace Foam::functionObjects::runTimeVis
{

Dictionaries::Dictionaries
    (
        const Time &runTime,
        const dictionary &postDict
    ) :
    hexMeshDict_(getHexMeshDict(runTime)),
    geometryDict_(getGeometryDict(runTime)),
    blockMeshDict_(getBlockMeshDict(runTime)),
    controlDict_(getControlDict(runTime)),
    meshObjectsDict_(getMeshObjectsDict(runTime)),
    rtppDict_(getRtppDict(runTime, postDict))
{}

const dictionary& Dictionaries::getRtppDict() const {return rtppDict_;}
const dictionary& Dictionaries::getHexMeshDict() const {return hexMeshDict_;}
const dictionary& Dictionaries::getBlockMeshDict() const {return blockMeshDict_;}
const dictionary& Dictionaries::getControlDict() const {return controlDict_;}
const dictionary& Dictionaries::getFunctionObjectsDict() const
{
    return controlDict_.subDict(controlDictKeys::FUNCTIONS_SUBDICT_KEY);
}
const dictionary& Dictionaries::getMeshObjectsDict() const {return meshObjectsDict_;}
const dictionary& Dictionaries::getGeometryDict() const {return geometryDict_;}

dictionary Dictionaries::getHexMeshDict(const Time &runTime)
{
    return getOptionalDictionary(runTime, HELYXHESHMESHDICT_NAME);
}

dictionary Dictionaries::getGeometryDict(const Time &runTime)
{
    return getOptionalDictionary(runTime, GEOMETRYDICT_NAME);
}

dictionary Dictionaries::getRtppDict(const Time &runTime, const dictionary& postDict)
{
    dictionary rtppDict = getOptionalDictionary(runTime, RTPP_DICT_NAME);
    if (rtppDict.empty())
    {
        return postDict;
    }
    else
    {
        return rtppDict;
    }
}

dictionary Dictionaries::getRequiredDictionary(const Time &runTime, const char* name)
{
    const fileName path = ParallelUtils::isRunningInParallel() ?
                          runTime.path() / runTime.caseSystem() / name :
                          runTime.systemPath() / name;
    return getDictionary(runTime, path);
}

dictionary Dictionaries::getDictionary(const Time &runTime, const fileName &path)
{
    return { /* NOLINT(cppcoreguidelines-slicing) */
        IOdictionary(
            IOobject(
                path,
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    };
}

dictionary Dictionaries::getBlockMeshDict(const Time &runTime)
{
    return getOptionalDictionary(runTime, BLOCKMESHDICT_NAME);
}

dictionary Dictionaries::getOptionalDictionary(const Time &runTime, const char* name)
{
    const fileName path = ParallelUtils::isRunningInParallel() ?
                          runTime.path() / runTime.caseSystem() / name :
                          runTime.systemPath() / name;

    if (IFstream(path).good())
    {
        return getDictionary(runTime, path);
    }
    else
    {
        return {};
    }
}

dictionary Dictionaries::getControlDict(const Time &runTime)
{
    return getRequiredDictionary(runTime, CONTROL_DICT_NAME);
}

dictionary Dictionaries::getMeshObjectsDict(const Time &runTime)
{
    return getOptionalDictionary(runTime, MESH_OBJECTS_DICT_NAME);
}

Dictionaries::~Dictionaries() = default;

}

// ************************************************************************* //
