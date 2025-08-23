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
    (c) 2013-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "linearLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::linearLaw> Foam::linearLaw::clone() const
{
    return autoPtr<linearLaw>(new linearLaw(*this));
}


inline Foam::autoPtr<Foam::linearLaw> Foam::linearLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<linearLaw>
    (
        new linearLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::linearLaw::rho(const scalar p) const
{
    return rho0_ + psi_*p;
}


inline Foam::scalar Foam::linearLaw::hContribution(scalar p) const
{
    // log(rho(p, T))/psi_
    return (p - Pstd)/rho(p);
}


inline Foam::scalar Foam::linearLaw::spContribution
(
    const scalar p,
    const scalar T
) const
{
    return -log((rho0_ + psi_*p)/(rho0_ + psi_*Pstd))/(T*psi_);
}


inline Foam::scalar Foam::linearLaw::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(rho(p)*R*T);
}


inline Foam::scalar Foam::linearLaw::psi() const
{
    return psi_;
}


// ************************************************************************* //
