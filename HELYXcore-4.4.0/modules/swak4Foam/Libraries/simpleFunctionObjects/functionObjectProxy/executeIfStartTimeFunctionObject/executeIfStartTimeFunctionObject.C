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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2011-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfStartTimeFunctionObject.H"
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
    defineTypeNameAndDebug(executeIfStartTimeFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfStartTimeFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfStartTimeFunctionObject::executeIfStartTimeFunctionObject
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
    ),
    runIfStartTime_(
        readBool(
            dict.lookup("runIfStartTime")
        )
    ),
    executeOnRestart_(
        dict.lookupOrDefault<bool>("executeOnRestart",false)
    ),
    timeIndexRead_(
        time().timeIndex()
    )
{
    Dbug<< " constructing " << name << endl;
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
    if (!dict.found("executeOnRestart")) {
        WarningIn("executeIfStartTimeFunctionObject::executeIfStartTimeFunctionObject")
            << "No entry 'executeOnRestart' in " << dict.name()
                << ". Assuming 'false'" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool executeIfStartTimeFunctionObject::condition()
{
    if (
        (time().timeIndex()==0 && !executeOnRestart_)
        ||
        (time().timeIndex()==timeIndexRead_ && executeOnRestart_)
    ) {
        return runIfStartTime_;
    } else {
        return !runIfStartTime_;
    }
}

bool executeIfStartTimeFunctionObject::read(const dictionary& dict)
{
    runIfStartTime_=
        readBool(
            dict.lookup("runIfStartTime")
        );

    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
