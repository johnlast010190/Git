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
    (c) 2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "readFields/readFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(readFields, 0);
    addToRunTimeSelectionTable(functionObject, readFields, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::readFields::readFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(),
    readOnStart_(true)
{
    read(dict);

    if (readOnStart_)
    {
        execute();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::readFields::~readFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::readFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    fieldSet_ = dict.lookup<wordList>("fields");
    dict.readIfPresent("readOnStart", readOnStart_);

    return true;
}


bool Foam::functionObjects::readFields::execute()
{
    forAll(fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];

        // If necessary load field
        loadField<scalar>(fieldName);
        loadField<vector>(fieldName);
        loadField<sphericalTensor>(fieldName);
        loadField<symmTensor>(fieldName);
        loadField<tensor>(fieldName);
    }

    return true;
}


bool Foam::functionObjects::readFields::write()
{
    return true;
}


// ************************************************************************* //
