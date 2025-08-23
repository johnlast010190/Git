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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Global
    readTimeControls

Description
    Read the control parameters used by setDeltaT

\*---------------------------------------------------------------------------*/

adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

//Table type entry for maxCo:
if (runTime.controlDict().found("maxCo"))
{
    maxCoDataPtr = Function1<scalar>::New("maxCo", runTime.controlDict());
    maxCo = maxCoDataPtr->value(runTime.timeIndex());
    if (maxCoDataPtr->type() == Function1Types::Constant<scalar>::typeName)
    {
        maxCoDataPtr.clear();
    }
}
else
{
    maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
}

maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

// ************************************************************************* //
