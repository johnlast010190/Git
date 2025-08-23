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
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hybridConstPerfectGas.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hybridConstPerfectGas>
Foam::hybridConstPerfectGas::clone() const
{
    return autoPtr<hybridConstPerfectGas>(new hybridConstPerfectGas(*this));
}


inline Foam::autoPtr<Foam::hybridConstPerfectGas>
Foam::hybridConstPerfectGas::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hybridConstPerfectGas>
    (
        new hybridConstPerfectGas(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hybridConstPerfectGas::rho
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return max(rhoMin_, min(rhoMax_, blend()*(p/(R*T)) + (1 - blend())*rho_));
}


inline Foam::scalar Foam::hybridConstPerfectGas::spContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar tRho = rho(p, T, R);
    if (tRho > rhoMax_ || tRho < rhoMin_)
    {
        return 0.0;
    }
    return blend()*(-R*log(p/Pstd));
}


inline Foam::scalar Foam::hybridConstPerfectGas::psi
(
    const scalar T,
    const scalar R
) const
{
    return blend()*(1.0/(R*T));
}


inline Foam::scalar Foam::hybridConstPerfectGas::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return blend()*1 + (1 - blend())*p/(rho_*R*T);
}


inline Foam::scalar Foam::hybridConstPerfectGas::CpMCv(const scalar R) const
{
    return blend()*R;
}


inline Foam::scalar Foam::hybridConstPerfectGas::hContribution(scalar p) const
{
    return (blend() - 1)*((p - Pstd)/rho_);
}


inline Foam::scalar Foam::hybridConstPerfectGas::alphav(const scalar T) const
{
    return blend()/T;
}


// ************************************************************************* //
