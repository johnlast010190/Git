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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "fieldExpression/fieldExpression.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldExpression, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldExpression::setResultName
(
    const word& typeName,
    const word& defaultArg
)
{
    if (fieldName_.empty())
    {
        fieldName_ = defaultArg;
    }

    if (resultName_.empty())
    {
        if (fieldName_ != defaultArg)
        {
            resultName_ = typeName + '(' + fieldName_ + ')';
        }
        else
        {
            resultName_ = typeName;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExpression::fieldExpression
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& fieldName,
    const word& resultName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(fieldName),
    resultName_(resultName),
    valueStr_(string::null)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExpression::~fieldExpression()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldExpression::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (fieldName_.empty() || dict.found("field"))
    {
        fieldName_ = dict.lookup<word>("field");
    }

    if (dict.found("result"))
    {
        resultName_ = dict.lookup<word>("result");
    }

    if (dict.found("value"))
    {
        valueStr_ = dict.lookup<string>("value");
    }

    return true;
}


bool Foam::functionObjects::fieldExpression::execute()
{
    if (!calc())
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        // Clear the result field from the objectRegistry if present
        clear();

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::fieldExpression::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::fieldExpression::clear()
{
    return clearObject(resultName_);
}


// ************************************************************************* //
