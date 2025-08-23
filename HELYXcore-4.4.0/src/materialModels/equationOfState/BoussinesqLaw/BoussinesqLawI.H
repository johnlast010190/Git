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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "BoussinesqLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::BoussinesqLaw> Foam::BoussinesqLaw::clone() const
{
    return autoPtr<BoussinesqLaw>(new BoussinesqLaw(*this));
}


inline Foam::autoPtr<Foam::BoussinesqLaw> Foam::BoussinesqLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<BoussinesqLaw>
    (
        new BoussinesqLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::BoussinesqLaw::rho(const scalar T) const
{
    return rho0_*(1.0 - beta_*(T - T0_));
}


inline Foam::scalar Foam::BoussinesqLaw::hContribution(scalar p, scalar T) const
{
    return (p - Pstd)/rho(T);
}


inline Foam::scalar Foam::BoussinesqLaw::CpContribution(scalar p, scalar T) const
{
    return rho0_*beta_*(p - Pstd)/sqr(rho(T));
}


inline Foam::scalar Foam::BoussinesqLaw::Z(scalar p, scalar T, scalar R) const
{
    return p/(rho(T)*R*T);
}


inline Foam::scalar Foam::BoussinesqLaw::alphav(scalar T) const
{
    return rho0_/rho(T)*beta_;
}

// ************************************************************************* //
