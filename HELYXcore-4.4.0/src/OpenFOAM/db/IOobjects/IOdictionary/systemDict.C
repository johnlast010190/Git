/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOobjects/IOdictionary/systemDict.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::IOobject Foam::systemDictIO
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName
)
{
    fileName dictPath = dictName;

    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];

        if
        (
            isDir
            (
                dictPath.isAbsolute()
              ? dictPath
              : ob.time().globalPath()/dictPath
            )
        )
        {
            dictPath = dictPath/dictName;
        }
    }

    Info<< "Reading " << dictPath << nl << endl;

    if (args.optionFound("dict") && !dictPath.isName())
    {
        return
            IOobject
            (
                dictPath.isAbsolute()
              ? dictPath
              : ob.time().globalPath()/dictPath,
                ob,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );
    }
    else
    {
        return
            IOobject
            (
                dictPath,
                ob.time().system(),
                regionName == polyMesh::defaultRegion ? word::null : regionName,
                ob,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );
    }
}


Foam::IOdictionary Foam::systemDict
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName
)
{
    return IOdictionary(systemDictIO(dictName, args, ob, regionName));
}


// ************************************************************************* //
