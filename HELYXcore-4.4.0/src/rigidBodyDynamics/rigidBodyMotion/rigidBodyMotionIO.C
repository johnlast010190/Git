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

#include "rigidBodyMotion/rigidBodyMotion.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::RBD::rigidBodyMotion::read(const dictionary& dict)
{
    rigidBodyModel::read(dict);

    aRelax_ = dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0);
    aDamp_ = dict.lookupOrDefault<scalar>("accelerationDamping", 1.0);
    report_ = dict.lookupOrDefault<Switch>("report", false);

    return true;
}


void Foam::RBD::rigidBodyMotion::write(Ostream& os) const
{
    rigidBodyModel::write(os);

    os.writeEntry("accelerationRelaxation", aRelax_);
    os.writeEntry("accelerationDamping", aDamp_);
    os.writeEntry("report", report_);
}


// ************************************************************************* //
