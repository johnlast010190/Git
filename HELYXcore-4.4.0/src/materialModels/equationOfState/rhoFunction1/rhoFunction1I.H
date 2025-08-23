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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rhoFunction1.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::rhoFunction1> Foam::rhoFunction1::clone() const
{
    return autoPtr<rhoFunction1>(new rhoFunction1(*this));
}


inline Foam::autoPtr<Foam::rhoFunction1> Foam::rhoFunction1::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<rhoFunction1>
    (
        new rhoFunction1(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::rhoFunction1::rho(const scalar T) const
{
    return max(rhoMin_, min(rhoMax_, rhoFunction1_->value(T)));
}


inline Foam::scalar Foam::rhoFunction1::CpContribution(scalar p, scalar T) const
{
    // At limits rho is constant function hence derivitive should be zero
    const scalar rho = rhoFunction1_->value(T);
    if (rho > rhoMax_ || rho < rhoMin_)
    {
        return 0.0;
    }

    return -((p - Pstd)/sqr(rho))*rhoFunction1_->derivative(T);
}


inline Foam::scalar Foam::rhoFunction1::hContribution(scalar p, scalar T) const
{
    return (p - Pstd)/rho(T);
}


inline Foam::scalar Foam::rhoFunction1::Z(scalar p, scalar T, scalar R) const
{
    return p/(rho(T)*R*T);
}


inline Foam::scalar Foam::rhoFunction1::alphav(scalar T) const
{
    return rhoFunction1_->derivative(T)/rho(T);
}

// ************************************************************************* //
