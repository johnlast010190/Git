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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "fieldValues/fieldValue/fieldValue.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldValue, 0);
    defineRunTimeSelectionTable(fieldValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& valueType
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, valueType, dict),
    scaleFactor_(1.0),
    dict_(dict),
    regionName_(word::null)
{
    read(dict);
}


Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const word& valueType
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(obr_, name, valueType, dict),
    scaleFactor_(1.0),
    dict_(dict),
    regionName_(word::null)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValue::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;
    }

    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fields_ = dict.lookup<wordList>("fields");
    writeFields_ = dict.lookup<Switch>("writeFields");
    scaleFactor_ = dict.lookupOrDefault<scalar>("scaleFactor", 1.0);

    return true;
}


bool Foam::functionObjects::fieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::fieldValue::write()
{
    Log << type() << " " << name() << " write:" << nl;

    return true;
}


// ************************************************************************* //
