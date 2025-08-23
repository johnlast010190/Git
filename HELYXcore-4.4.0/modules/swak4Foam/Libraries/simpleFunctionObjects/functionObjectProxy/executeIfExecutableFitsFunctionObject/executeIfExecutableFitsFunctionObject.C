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
    2011-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfExecutableFitsFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"


#ifdef darwin
#include "mach-o/dyld.h"
#endif
#ifdef __linux__
#include <unistd.h>
#endif
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(executeIfExecutableFitsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfExecutableFitsFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfExecutableFitsFunctionObject::executeIfExecutableFitsFunctionObject
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
    readRegexp(dict);

    fileName exePath;

#ifdef darwin
    {
        char path[1024];
        uint32_t size = sizeof(path);
        if (_NSGetExecutablePath(path, &size) == 0) {
            exePath=string(path);
        }
    }
#elif defined(__linux__)
    {
        const int bufSize=1024;
        char path[bufSize];
        label length=readlink("/proc/self/exe",path,bufSize-1);
        path[length]='\0';
        exePath=string(path);
    }
#else
    NotImplemented;
#endif

    executable_=exePath.name();

    if (debug) {
        Info<< "Executable: " << executable_ << " "<< exePath << endl;
    }

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool executeIfExecutableFitsFunctionObject::condition()
{
    return executableNameRegexp_.match(executable_);
}

void executeIfExecutableFitsFunctionObject::readRegexp(const dictionary& dict)
{
    executableNameRegexp_.set(
        string(dict.lookup("executableNameRegexp")),
        dict.lookupOrDefault<bool>("ignoreCase",false)
    );
}

bool executeIfExecutableFitsFunctionObject::read(const dictionary& dict)
{
    readRegexp(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
