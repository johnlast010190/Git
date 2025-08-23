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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline Foam::scalar Foam::TimeState::timeOutputValue() const
{
    return timeToUserTime(value());
}


inline Foam::label Foam::TimeState::timeIndex() const
{
    return timeIndex_;
}


inline Foam::scalar Foam::TimeState::deltaTValue() const
{
    return deltaT_;
}


inline Foam::scalar Foam::TimeState::deltaT0Value() const
{
    return deltaT0_;
}


inline Foam::dimensionedScalar Foam::TimeState::deltaT() const
{
    return dimensionedScalar("deltaT", dimTime, deltaT_);
}


inline Foam::dimensionedScalar Foam::TimeState::deltaT0() const
{
    return dimensionedScalar("deltaT0", dimTime, deltaT0_);
}


inline bool Foam::TimeState::writeTime() const
{
    return writeTime_;
}


inline bool Foam::TimeState::outputTime() const
{
    return writeTime_;
}


// ************************************************************************* //
