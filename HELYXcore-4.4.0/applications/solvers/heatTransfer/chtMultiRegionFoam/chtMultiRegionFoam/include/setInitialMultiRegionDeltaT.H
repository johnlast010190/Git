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
    (c) 2020 Engys Ltd.

Global
    setInitialMultiRegionDeltaT

Description
    Set the initial timestep for the CHT MultiRegion solver.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && ((CoNum > SMALL) || (DiNum > SMALL)))
    {
        if (CoNum < SMALL)
        {
            CoNum = SMALL;
        }

        if (DiNum < SMALL)
        {
            DiNum = SMALL;
        }

        runTime.setDeltaT
        (
            min
            (
                min(maxCo/CoNum, maxDi/DiNum)*runTime.deltaT().value(),
                min(runTime.deltaTValue(), maxDeltaT)
            )
        );
        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;

        //exit Helyx if deltaT is below user secified minDeltaT
        if (runTime.deltaTValue() < minDeltaT)
        {
            FatalErrorInFunction
                << "minDeltaT = " << minDeltaT << nl
                << "The computed deltaT is lower than the minDeltaT"
                << exit(FatalError);
        }

    }
}

// ************************************************************************* //
