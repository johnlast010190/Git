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
    2012-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "writeAndEndFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeAndEndFunctionObject, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeAndEndFunctionObject::writeAndEndFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    isStopped_(false),
    storeAndWritePreviousState_(
        dict.lookupOrDefault<bool>("storeAndWritePreviousState",false)
    )
{
    if (storeAndWritePreviousState_) {
        lastTimes_.set(
            new TimeCloneList(
                dict
            )
        );
    } else if (!dict.found("storeAndWritePreviousState")){
        WarningIn("writeAndEndFunctionObject::writeAndEndFunctionObject")
            << "'storeAndWritePreviousState' unset in " << dict.name()
                << endl;
    }
}

bool writeAndEndFunctionObject::start()
{
    if (debug) {
        Info<< name() << "::start() - Entering" << endl;
    }

    simpleFunctionObject::start();

    if (debug) {
        Info<< name() << "::start() - Leaving" << endl;
    }

    return true;
}

void writeAndEndFunctionObject::writeSimple()
{
    if (debug) {
        Info<< name() << "::writeSimple() - Entering" << endl;
    }
    if (isStopped()) {
        if (debug) {
            Info<< name() << "::writeSimple() - isStopped" << endl;
        }
        return;
    }
    if (
        this->endRunNow()
    ) {
        if (debug) {
            Info<< name() << "::writeSimple() - stopping" << endl;
        }
        isStopped_=true;

        Info<< "Ending run because of functionObject " << this->name() << endl;
        Info<< "Writing current state" << endl;
        const_cast<Time &>(time()).writeAndEnd();
        if (lastTimes_.valid()) {
            Pout<< "Writing old times" << endl;
            lastTimes_->write();
        }
    }
    if (lastTimes_.valid()) {
        Dbug<< "Storing current time" << endl;
        lastTimes_->copy(time());
    }
    if (debug) {
        Info<< name() << "::writeSimple() - Leaving" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
