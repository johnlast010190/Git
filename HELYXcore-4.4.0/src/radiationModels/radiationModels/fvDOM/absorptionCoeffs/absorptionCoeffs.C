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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationModels/fvDOM/absorptionCoeffs/absorptionCoeffs.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionCoeffs::~absorptionCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::absorptionCoeffs::checkT(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "using absorptionCoeffs out of temperature range:" << nl
            << "    " << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << nl << endl;
    }
}


const Foam::radiationModels::absorptionCoeffs::coeffArray&
Foam::radiationModels::absorptionCoeffs::coeffs
(
    const scalar T
) const
{
    checkT(T);

    if (T < Tcommon_)
    {
        return lowACoeffs_;
    }
    else
    {
        return highACoeffs_;
    }
}


void Foam::radiationModels::absorptionCoeffs::initialise(const dictionary& dict)
{
    Tcommon_ = dict.lookup<scalar>("Tcommon");
    Tlow_ = dict.lookup<scalar>("Tlow");
    Thigh_ = dict.lookup<scalar>("Thigh");
    invTemp_ = dict.lookup<bool>("invTemp");
    lowACoeffs_ = dict.lookup<coeffArray>("loTcoeffs");
    highACoeffs_ = dict.lookup<coeffArray>("hiTcoeffs");
}


// ************************************************************************* //
