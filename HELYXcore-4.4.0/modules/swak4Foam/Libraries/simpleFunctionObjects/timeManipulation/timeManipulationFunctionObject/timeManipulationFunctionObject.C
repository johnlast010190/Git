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
    2008-2011, 2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "timeManipulationFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeManipulationFunctionObject, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeManipulationFunctionObject::timeManipulationFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    tolerateAdaptiveTimestep_(false),
    myEndTime_(-1)
{
}

bool timeManipulationFunctionObject::start()
{
    simpleFunctionObject::start();

    tolerateAdaptiveTimestep_=dict_.lookupOrDefault<bool>(
        "tolerateAdaptiveTimestep",false
    );

    if (
        !tolerateAdaptiveTimestep_
        &&
        time().controlDict().lookupOrDefault<bool>("adjustTimeStep",false)
    ) {
        FatalErrorIn("timeManipulationFunctionObject::start()")
            << "'adjustTimeStep' is set. Function object " << name()
                << " manipulates time. This may lead to strange behaviour. "
                << " If you think that is OK set 'tolerateAdaptiveTimestep'"
                << " in " << name()
                << endl
                << exit(FatalError);
    }

    return true;
}

void timeManipulationFunctionObject::writeSimple()
{
    scalar newDeltaT=this->deltaT();
    scalar newEndTime=this->endTime();

    scalar minDeltaT=newDeltaT;
    scalar maxDeltaT=newDeltaT;
    reduce(
        std::tie(minDeltaT, maxDeltaT),
        ParallelOp<minOp<scalar>, maxOp<scalar>>{}
    );
    if (minDeltaT!=maxDeltaT) {
        FatalErrorIn("timeManipulationFunctionObject::writeSimple()")
            << "Across the processors the minimum " << minDeltaT
                << " and the maximum " << maxDeltaT << " of the new deltaT"
                << " differ by " << maxDeltaT-minDeltaT
                << endl
                << exit(FatalError);

    }
    if (newDeltaT!=time().deltaT().value()) {
        Info<< name() << " sets deltaT from " <<
            time().deltaT().value() << " to " << newDeltaT << endl;

        const_cast<Time&>(time()).setDeltaT(newDeltaT);
    }

    scalar minEndTime=newEndTime;
    scalar maxEndTime=newEndTime;
    reduce(
        std::tie(minEndTime, maxEndTime),
        ParallelOp<minOp<scalar>, maxOp<scalar>>{}
    );
    if (minEndTime!=maxEndTime) {
        FatalErrorIn("timeManipulationFunctionObject::writeSimple()")
            << "Across the processors the minimum " << minEndTime
                << " and the maximum " << maxEndTime << " of the new endTime"
                << " differ by " << maxEndTime-minEndTime
                << endl
                << exit(FatalError);

    }
    if (newEndTime!=time().endTime().value()) {
        Info<< name() << " sets endTime from " <<
            time().endTime().value() << " to " << newEndTime << endl;

        const_cast<Time&>(time()).setEndTime(newEndTime);

        myEndTime_=newEndTime;
    }

    if (
        myEndTime_>-1
        &&
        time().value()>=time().endTime().value()
        &&
        !time().outputTime()
    ) {
        WarningIn("timeManipulationFunctionObject::writeSimple()")
            << "Forcing write because we (" << name()
                << ") changed the endTime to "
                << myEndTime_ << " and this is not a write-time"
                << endl;

        const_cast<Time&>(time()).writeNow();
    }
}

scalar timeManipulationFunctionObject::deltaT()
{
    return time().deltaT().value();
}

scalar timeManipulationFunctionObject::endTime()
{
    return time().endTime().value();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
