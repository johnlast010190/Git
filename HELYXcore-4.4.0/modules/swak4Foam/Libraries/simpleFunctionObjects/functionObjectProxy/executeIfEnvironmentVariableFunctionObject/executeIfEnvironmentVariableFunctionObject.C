/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfEnvironmentVariableFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"
#include "include/OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(executeIfEnvironmentVariableFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfEnvironmentVariableFunctionObject,
        dictionary
    );

    template<>
    const char* NamedEnum<Foam::executeIfEnvironmentVariableFunctionObject::fitVariableMode,3>::names[]=
    {
        "exists",
        "doesNotExist",
        "fitsRegexp"
    };
    const NamedEnum<executeIfEnvironmentVariableFunctionObject::fitVariableMode,3> executeIfEnvironmentVariableFunctionObject::fitVariableModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfEnvironmentVariableFunctionObject::executeIfEnvironmentVariableFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    conditionalFunctionObjectListProxy(
        name,
        t,
        dict
    )
{
    // do it here to avoid the superclass-read being read twice
    readData(dict);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool executeIfEnvironmentVariableFunctionObject::condition()
{
    bool exists=env(variableName_);
    if (writeDebug()) {
        Info<< "Variable " << variableName_ << " "
            << (exists ? "exists" : "does not exist")
            << endl;
    }

    string content;

    switch(fitVariableMode_) {
        case fvmExists:
            return exists;
        case fvmDoesNotExist:
            return !exists;
        case fvmFitsRegexp:
            if (exists) {
                content=getEnv(variableName_);
                if (writeDebug()) {
                    Info<< "Content of " << variableName_ << ": "
                        << content
                        << (contentRegexp_.match(content) ? " fits " : " does not fit ")
                        << endl;
                }
                return contentRegexp_.match(content);
            }
            break;
        default:
            FatalErrorIn("executeIfEnvironmentVariableFunctionObject::condition()")
                << "The fitting mode " << fitVariableModeNames_[fitVariableMode_]
                    << " is not implemented" << endl
                    << exit(FatalError);
    }
    return false;
}

void executeIfEnvironmentVariableFunctionObject::readData(const dictionary& dict)
{
    fitVariableMode_=fitVariableModeNames_[dict.lookup("fitVariableMode")];
    variableName_=word(dict.lookup("variableName"));

    if (fitVariableMode_==fvmFitsRegexp) {
        contentRegexp_.set(
            string(dict.lookup("contentRegexp")),
            dict.lookupOrDefault<bool>("ignoreCase",false)
        );
    }
}

bool executeIfEnvironmentVariableFunctionObject::read(const dictionary& dict)
{
    readData(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
