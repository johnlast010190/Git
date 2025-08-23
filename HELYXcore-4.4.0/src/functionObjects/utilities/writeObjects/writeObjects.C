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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "writeObjects/writeObjects.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeObjects, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeObjects,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::writeObjects::writeOption,
    3
>::names[] =
{
    "autoWrite",
    "noWrite",
    "anyWrite"
};

const Foam::NamedEnum
<
    Foam::functionObjects::writeObjects::writeOption,
    3
> Foam::functionObjects::writeObjects::writeOptionNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::writeObjects
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    writeOption_(ANY_WRITE),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::~writeObjects()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeObjects::read(const dictionary& dict)
{
    functionObject::read(dict);

    if (dict.found("field"))
    {
        objectNames_.setSize(1);
        objectNames_[0] = dict.lookup<wordRe>("field");
    }
    else if (dict.found("fields"))
    {
        objectNames_ = dict.lookup<wordReList>("fields");
    }
    else
    {
        objectNames_ = dict.lookup<wordReList>("objects");
    }

    if (dict.found("writeOption"))
    {
        writeOption_ = writeOptionNames_.read(dict.lookup("writeOption"));
    }
    else
    {
        writeOption_ = ANY_WRITE;
    }

    return true;
}


bool Foam::functionObjects::writeObjects::execute()
{
    return true;
}


bool Foam::functionObjects::writeObjects::write()
{
    Log << type() << " " << name() << " write:" << nl;

    if (!obr_.time().writeTime())
    {
        obr_.time().writeTimeDict();
    }

    DynamicList<word> allNames(obr_.toc().size());
    forAll(objectNames_, i)
    {
        wordList names(obr_.names<regIOobject>(objectNames_[i]));

        if (names.size())
        {
            allNames.append(names);
        }
        else
        {
            WarningInFunction
                << "Object " << objectNames_[i] << " not found in "
                << "database. Available objects:" << nl << obr_.sortedToc()
                << endl;
        }
    }

    forAll(allNames, i)
    {
        regIOobject& obj = const_cast<regIOobject&>
        (
            obr_.lookupObject<regIOobject>(allNames[i])
        );

        switch (writeOption_)
        {
            case AUTO_WRITE:
            {
                if (obj.writeOpt() != IOobject::AUTO_WRITE)
                {
                    continue;
                }

                break;
            }
            case NO_WRITE:
            {
                if (obj.writeOpt() != IOobject::NO_WRITE)
                {
                    continue;
                }

                break;
            }
            case ANY_WRITE:
            {
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown writeOption "
                    << writeOptionNames_[writeOption_]
                    << ". Valid writeOption types are "
                    << writeOptionNames_
                    << exit(FatalError);
            }
        }

        if
        (
            obj.writeOpt() == IOobject::AUTO_WRITE
         && obr_.time().writeTime()
        )
        {
            Log << "    automatically written object " << obj.name() << endl;
        }
        else
        {
            Log << "    writing object " << obj.name() << endl;

            obj.write();
        }
    }

    return true;
}


// ************************************************************************* //
