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
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/RHWindProfile/RHWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(RHWindKineticEnergyProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindKineticEnergyProfile::RHWindKineticEnergyProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    RHWindProfile(entryName,dict)
{
}


Foam::Function1Types::RHWindKineticEnergyProfile::RHWindKineticEnergyProfile
(
     const RHWindKineticEnergyProfile& wp
)
:
    RHWindProfile(wp)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindKineticEnergyProfile::~RHWindKineticEnergyProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::RHWindKineticEnergyProfile::value(const scalar z_) const
{
    return magSqr(Foam::Function1Types::RHWindProfile::Ufrictional())/sqrt(Cmu_);
}
// ************************************************************************* //
