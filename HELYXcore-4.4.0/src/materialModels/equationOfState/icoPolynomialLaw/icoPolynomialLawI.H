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
    (c) 2016-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "icoPolynomialLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::icoPolynomialLaw>
Foam::icoPolynomialLaw::clone() const
{
    return autoPtr<icoPolynomialLaw>(new icoPolynomialLaw(*this));
}


inline Foam::autoPtr<Foam::icoPolynomialLaw> Foam::icoPolynomialLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<icoPolynomialLaw>
    (
        new icoPolynomialLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::icoPolynomialLaw::rho(const scalar T) const
{
    return max(rhoMin_, min(rhoMax_, rhoCoeffs_.value(T)));
}


inline Foam::scalar Foam::icoPolynomialLaw::hContribution(scalar p, scalar T) const
{
    return (p - Pstd)/rho(T);
}


inline Foam::scalar Foam::icoPolynomialLaw::CpContribution(scalar p, scalar T) const
{
    // At limits rho is constant function hence derivitive should be zero
    const scalar rho = rhoCoeffs_.value(T);
    if (rho > rhoMax_ || rho < rhoMin_)
    {
        return 0.0;
    }

    return -((p - Pstd)/sqr(rho))*rhoCoeffs_.derivative(T);
}


inline Foam::scalar Foam::icoPolynomialLaw::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(rho(T)*R*T);
}


inline Foam::scalar Foam::icoPolynomialLaw::alphav(const scalar T) const
{
    return -rhoCoeffs_.derivative(T)/rho(T);
}


// ************************************************************************* //
