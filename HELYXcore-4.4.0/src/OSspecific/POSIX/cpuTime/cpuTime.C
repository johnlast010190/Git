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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cpuTime/cpuTime.H"
#include <unistd.h>

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const long Foam::cpuTime::Hz_(sysconf(_SC_CLK_TCK));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cpuTime::getTime(timeType& t)
{
    times(&t);
}


double Foam::cpuTime::timeDifference(const timeType& beg, const timeType& end)
{
    return
    (
        double
        (
            (end.tms_utime + end.tms_stime)
          - (beg.tms_utime + beg.tms_stime)
        )/Hz_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuTime::cpuTime()
{
    getTime(startTime_);
    lastTime_ = startTime_;
    newTime_ = startTime_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::cpuTime::elapsedCpuTime() const
{
    getTime(newTime_);
    return timeDifference(startTime_, newTime_);
}


double Foam::cpuTime::cpuTimeIncrement() const
{
    lastTime_ = newTime_;
    getTime(newTime_);
    return timeDifference(lastTime_, newTime_);
}


// ************************************************************************* //
