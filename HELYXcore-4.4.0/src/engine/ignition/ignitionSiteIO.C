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

#include "ignition/ignitionSite.H"
#include "engineTime/engineTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ignitionSite::ignitionSite
(
    Istream& is,
    const Time& db,
    const fvMesh& mesh
)
:
    db_(db),
    mesh_(mesh),
    ignitionSiteDict_(is),
    location_(ignitionSiteDict_.lookup("location")),
    diameter_(ignitionSiteDict_.lookup<scalar>("diameter")),
    time_
    (
        db_.userTimeToTime
        (
            ignitionSiteDict_.lookup<scalar>("start")
        )
    ),
    duration_
    (
        db_.userTimeToTime
        (
            ignitionSiteDict_.lookup<scalar>("duration")
        )
    ),
    strength_(ignitionSiteDict_.lookup<scalar>("strength")),
    timeIndex_(db_.timeIndex())
{
    // Check state of Istream
    is.check(FUNCTION_NAME);

    findIgnitionCells(mesh_);
}


Foam::ignitionSite::ignitionSite
(
    Istream& is,
    const engineTime& edb,
    const fvMesh& mesh
)
:
    db_(edb),
    mesh_(mesh),
    ignitionSiteDict_(is),
    location_(ignitionSiteDict_.lookup("location")),
    diameter_(ignitionSiteDict_.lookup<scalar>("diameter")),
    time_
    (
        db_.userTimeToTime
        (
            edb.degToTime(ignitionSiteDict_.lookup<scalar>("start"))
        )
    ),
    duration_
    (
        db_.userTimeToTime
        (
            edb.degToTime(ignitionSiteDict_.lookup<scalar>("duration"))
        )
    ),
    strength_(ignitionSiteDict_.lookup<scalar>("strength")),
    timeIndex_(db_.timeIndex())
{
    // Check state of Istream
    is.check(FUNCTION_NAME);

    findIgnitionCells(mesh_);
}


// ************************************************************************* //
