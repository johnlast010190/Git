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
    (c) 2019-2024 Engys Ltd.
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fieldsExpression/fieldsExpression.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldsExpression, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldsExpression::setResultName
(
    const word& typeName,
    const wordList& defaultArgs
)
{
    if (fieldNames_.empty())
    {
        fieldNames_ = defaultArgs;
    }

    if (resultName_.empty())
    {
        if (!fieldNames_.empty())
        {
            resultName_ = typeName + '(' + fieldNames_[0];
            for (label i=1; i<fieldNames_.size(); i++)
            {
                resultName_ += ',' + fieldNames_[i];
            }

            if (!valueStr_.empty())
            {
                resultName_ += ',' + name();
            }

            resultName_ += ')';
        }
        else
        {
            resultName_ = typeName;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldsExpression::fieldsExpression
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const wordList& fieldNames,
    const word& resultName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(fieldNames),
    resultName_(resultName),
    valueStr_(string::null)
{
    read(dict);

    if (fieldNames_.size() < 2 && valueStr_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "functionObject::" << type() << " " << name
            << " requires at least 2 fields (if value entry not set). Only "
            << fieldNames_.size() << " provided: " << fieldNames_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldsExpression::~fieldsExpression()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldsExpression::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (fieldNames_.empty() || dict.found("fields"))
    {
        fieldNames_ = dict.lookup<wordList>("fields");
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


bool Foam::functionObjects::fieldsExpression::execute()
{
    if (!calc())
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " cannot find required fields " << fieldNames_ << endl;

        // Clear the result fields from the objectRegistry if present
        clear();

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::fieldsExpression::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::fieldsExpression::clear()
{
    return clearObject(resultName_);
}


// ************************************************************************* //
