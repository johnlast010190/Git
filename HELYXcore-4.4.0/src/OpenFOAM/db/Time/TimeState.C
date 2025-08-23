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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimeState::TimeState()
:
    dimensionedScalar(Time::timeName(0), dimTime, 0),
    timeIndex_(0),
    deltaT_(0),
    deltaTDataPtr_(),
    deltaTSave_(0),
    deltaT0_(0),
    deltaTchanged_(false),
    writeTimeIndex_(0),
    writeTime_(false),
    deltaTcontrol_()
{}


Foam::TimeState::TimeState(const TimeState& ts)
:
    dimensionedScalar(ts),
    timeIndex_(ts.timeIndex_),
    deltaT_(ts.deltaT_),
    deltaTDataPtr_(ts.deltaTDataPtr_, false),
    deltaTSave_(ts.deltaTSave_),
    deltaT0_(ts.deltaT0_),
    deltaTchanged_(ts.deltaTchanged_),
    writeTimeIndex_(ts.writeTimeIndex_),
    writeTime_(ts.writeTime_),
    deltaTcontrol_(ts.deltaTcontrol_)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::TimeState::~TimeState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::TimeState::userTimeToTime(const scalar theta) const
{
    return theta;
}


Foam::scalar Foam::TimeState::timeToUserTime(const scalar t) const
{
    return t;
}


Foam::TimeState& Foam::TimeState::operator=(const TimeState& ts)
{
    dimensionedScalar::operator=(ts);
    timeIndex_ = ts.timeIndex_;
    deltaT_ = ts.deltaT_;
    if (ts.deltaTDataPtr_.valid())
    {
        deltaTDataPtr_.reset(ts.deltaTDataPtr_().clone().ptr());
    }
    else
    {
        deltaTDataPtr_.clear();
    }
    deltaTSave_ = ts.deltaTSave_;
    deltaT0_ = ts.deltaT0_;
    deltaTchanged_ = ts.deltaTchanged_;
    writeTimeIndex_ = ts.writeTimeIndex_;
    writeTime_ = ts.writeTime_;
    deltaTcontrol_ = ts.deltaTcontrol_;

    return *this;
}

// ************************************************************************* //
