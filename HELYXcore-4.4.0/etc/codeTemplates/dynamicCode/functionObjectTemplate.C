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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "functionObjectTemplate.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/unitConversion/unitConversion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(${typeName}FunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    ${typeName}FunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& ${typeName}FunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FunctionObject::${typeName}FunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FunctionObject::~${typeName}FunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ${typeName}FunctionObject::read(const dictionary& dict)
{
    if (${verbose:-false})
    {
        Info<<"read ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeRead}
//}}} end code

    return true;
}


bool ${typeName}FunctionObject::execute()
{
    if (${verbose:-false})
    {
        Info<<"execute ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeExecute}
//}}} end code

    return true;
}


bool ${typeName}FunctionObject::write()
{
    if (${verbose:-false})
    {
        Info<<"write ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeWrite}
//}}} end code

    return true;
}


bool ${typeName}FunctionObject::end()
{
    if (${verbose:-false})
    {
        Info<<"end ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeEnd}
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
