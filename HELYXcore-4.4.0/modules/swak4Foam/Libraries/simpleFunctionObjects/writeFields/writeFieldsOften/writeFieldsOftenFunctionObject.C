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

#include "writeFieldsOftenFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeFieldsOftenFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeFieldsOftenFunctionObject,
        dictionary
    );


    // copied from Time.C because the original is protected
    // to work the order of values in writeControls must not change
template<>
const char* NamedEnum<Foam::Time::writeControls, 6>::names[] =
{
    "timeStep",
    "runTime",
    "adjustableRunTime",
    "clockTime",
    "cpuTime",
    "baselineTime"
};

const NamedEnum<Foam::Time::writeControls, 6> writeControlNames;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeFieldsOftenFunctionObject::writeFieldsOftenFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    writeFieldsGeneralFunctionObject(name,t,dict),
    writeControl_(Time::wcTimeStep),
    writeInterval_(GREAT),
    outputTimeIndex_(0),
    baselineTime_(0)
{
}

bool writeFieldsOftenFunctionObject::start()
{
    writeFieldsGeneralFunctionObject::start();

    word wcName=dict_.lookup("writeControl");

    writeControl_ = writeControlNames.read
        (
            dict_.lookup("writeControl")
        );

    writeInterval_ = dict_.lookup<scalar>("writeIntervall");
    if (writeControl_ == Time::wcTimeStep && label(writeInterval_) <1) {
        WarningIn("bool writeFieldsOftenFunctionObject::start()")
            << "writeInterval " << writeInterval_
                << " < 1 for writeControl timeStep. Reseting to 1 "<< endl;
        writeInterval_=1;
    }
    baselineTime_ = dict_.lookup<scalar>("baselineTime");

    Info<< "Additional fields " << fieldNames() << " will be written "
        << "with writeControl " << wcName << " and intervall " << writeInterval_ << endl;

    if (writeControl_ == Time::wcAdjustableRunTime) {
        WarningIn("bool writeFieldsOftenFunctionObject::start()")
            << "Cant adjust the run-time. Defaulting to runTime" << endl;

    }

    outputTimeIndex_=0;

    return true;
}

bool writeFieldsOftenFunctionObject::outputTime(const bool forceWrite)
{
    if (forceWrite) {
        return true;
    }

    if (time().time().value()<after()) {
        return false;
    }

    bool writeNow=false;

    // lifted from Time::operator++
    switch(writeControl_)
    {
        case Time::wcTimeStep:
            writeNow = !(time().timeIndex()%label(writeInterval_));
            break;

        case Time::wcRunTime:
        case Time::wcAdjustableRunTime:
            {
                label outputTimeIndex =
                    label(((time().time().value() - time().startTime().value()) + 0.5*time().deltaT().value())/writeInterval_);

                if (outputTimeIndex > outputTimeIndex_)
                {
                    writeNow = true;
                    outputTimeIndex_ = outputTimeIndex;
                }
                else
                {
                    writeNow = false;
                }
            }
        break;

        case Time::wcCpuTime:
            {
                label outputTimeIndex =
                    label(time().elapsedCpuTime()/writeInterval_);

                if (outputTimeIndex > outputTimeIndex_)
                {
                    writeNow = true;
                    outputTimeIndex_ = outputTimeIndex;
                }
                else
                {
                    writeNow = false;
                }
            }
        break;

        case Time::wcClockTime:
            {
                label outputTimeIndex = label(time().elapsedClockTime()/writeInterval_);
                if (outputTimeIndex > outputTimeIndex_)
                {
                    writeNow = true;
                    outputTimeIndex_ = outputTimeIndex;
                }
                else
                {
                    writeNow = false;
                }
            }
        break;

        case Time::wcBaselineTime:
            {
                scalar tStart = time().time().value() - time().deltaT().value();

                if (time().time().value() >= baselineTime_ && tStart >= baselineTime_)
                {
                    label outputIndexStart = label
                    (
                        ((tStart - baselineTime_))
                        / writeInterval_
                    );

                    label outputIndexEnd = label
                    (
                        (time().time().value() - baselineTime_)
                        / writeInterval_
                    );

                    if (outputIndexEnd > outputIndexStart)
                    {
                        writeNow = true;
                        outputTimeIndex_ = outputIndexEnd;
                    }
                }
            }
        break;

        case Time::wcOptimizationCycle:
        {
            int intValue = round(time().time().value());
            const scalar timeTol =
                max
                (
                    min
                    (
                        pow(10.0, -time().time().getTimeNamePrecision()),
                        0.1*time().deltaT().value()
                    ),
                    SMALL
                );

            if (mag(time().time().value() - intValue) < timeTol)
            {
                writeNow = !(intValue % int (writeInterval_));
            }
        }
        break;
    };

    return writeNow;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
