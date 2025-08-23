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
    (c) 2019-2023 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rPolynomialLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::rPolynomialLaw>
Foam::rPolynomialLaw::clone() const
{
    return autoPtr<rPolynomialLaw>(new rPolynomialLaw(*this));
}


inline Foam::autoPtr<Foam::rPolynomialLaw> Foam::rPolynomialLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<rPolynomialLaw>
    (
        new rPolynomialLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::rPolynomialLaw::rho
(
    const scalar p,
    const scalar T
) const
{
    return 1/(C_[0] + (C_[1] + C_[2]*T - C_[4]*p)*T - C_[3]*p);
}


inline Foam::scalar Foam::rPolynomialLaw::hContribution
(
    scalar pAbs,
    scalar T
) const
{
    // hdep = p/rho - /int0^T (cp - cv) dT =
    //     p/rho - \int0^T((T*(C_[1] + 2C_[2]T - C_[4]p)^2)/(C_[3] + C_[4]T))dT
    // A lot of coefficients can be pre-computed on construction to
    // increase the performance of this function.
    const scalar p = pAbs - Pstd;
    return
        p/rho(p, T)
      - (
            4*sqr(C_[2])*pow(C_[4], 3.0)*pow(T, 3)/3
          + 2*C_[2]*sqr(C_[4])*sqr(T)*(C_[1]*C_[4] - C_[2]*C_[3] - sqr(C_[4])*p)
          + C_[3]*(log(C_[3]) - log(C_[3] + C_[4]*T))
           *sqr(-C_[1]*C_[4] + 2*C_[2]*C_[3] + sqr(C_[4])*p)
          + C_[4]*T
           *(
                sqr(C_[1])*sqr(C_[4])
              - 4*C_[1]*C_[2]*C_[3]*C_[4]
              - 2*C_[1]*pow(C_[4], 3)*p
              + 4*sqr(C_[2])*sqr(C_[3])
              + 4*C_[2]*C_[3]*sqr(C_[4])*p
              + pow(C_[4], 4)*sqr(p)
            )
        )/pow(C_[4], 4);
}


inline Foam::scalar Foam::rPolynomialLaw::eContribution
(
    scalar pAbs,
    scalar T
) const
{
    // Edep = Hdep - p/rho
    const scalar p = pAbs - Pstd;
    return
      - (
            4*sqr(C_[2])*pow(C_[4], 3.0)*pow(T, 3)/3
          + 2*C_[2]*sqr(C_[4])*sqr(T)*(C_[1]*C_[4] - C_[2]*C_[3] - sqr(C_[4])*p)
          + C_[3]*(log(C_[3]) - log(C_[3] + C_[4]*T))
           *sqr(-C_[1]*C_[4] + 2*C_[2]*C_[3] + sqr(C_[4])*p)
          + C_[4]*T
           *(
                sqr(C_[1])*sqr(C_[4])
              - 4*C_[1]*C_[2]*C_[3]*C_[4]
              - 2*C_[1]*pow(C_[4], 3)*p
              + 4*sqr(C_[2])*sqr(C_[3])
              + 4*C_[2]*C_[3]*sqr(C_[4])*p
              + pow(C_[4], 4)*sqr(p)
            )
        )/pow(C_[4], 4);
}


inline Foam::scalar Foam::rPolynomialLaw::CpContribution
(
    scalar pAbs,
    scalar T
) const
{
    // CpDep = dHdep/dT
    // A lot of coefficients can be pre-computed on construction to
    // increase the performance of this function.
    const scalar p = pAbs - Pstd;
    return
        p*(C_[1] + 2*C_[2]*T - C_[4]*p)
      - (
            12*sqr(C_[2])*pow(C_[4], 3.0)*pow(T, 2)/3
          + 4*C_[2]*sqr(C_[4])*T*(C_[1]*C_[4] - C_[2]*C_[3] - sqr(C_[4])*p)
          + -C_[3]*(C_[4]/(C_[3] + C_[4]*T))
           *sqr(-C_[1]*C_[4] + 2*C_[2]*C_[3] + sqr(C_[4])*p)
          + C_[4]
           *(
                sqr(C_[1])*sqr(C_[4])
              - 4*C_[1]*C_[2]*C_[3]*C_[4]
              - 2*C_[1]*pow(C_[4], 3)*p
              + 4*sqr(C_[2])*sqr(C_[3])
              + 4*C_[2]*C_[3]*sqr(C_[4])*p
              + pow(C_[4], 4)*sqr(p)
            )
        )/pow(C_[4], 4);
}


inline Foam::scalar Foam::rPolynomialLaw::CvContribution
(
    scalar pAbs,
    scalar T
) const
{
    // CvDep = dEdep/dT
    const scalar p = pAbs - Pstd;
    return
      - (
            12*sqr(C_[2])*pow(C_[4], 3.0)*pow(T, 2)/3
          + 4*C_[2]*sqr(C_[4])*T*(C_[1]*C_[4] - C_[2]*C_[3] - sqr(C_[4])*p)
          + -C_[3]*(C_[4]/(C_[3] + C_[4]*T))
           *sqr(-C_[1]*C_[4] + 2*C_[2]*C_[3] + sqr(C_[4])*p)
          + C_[4]
           *(
                sqr(C_[1])*sqr(C_[4])
              - 4*C_[1]*C_[2]*C_[3]*C_[4]
              - 2*C_[1]*pow(C_[4], 3)*p
              + 4*sqr(C_[2])*sqr(C_[3])
              + 4*C_[2]*C_[3]*sqr(C_[4])*p
              + pow(C_[4], 4)*sqr(p)
            )
        )/pow(C_[4], 4);
}


inline Foam::scalar Foam::rPolynomialLaw::Z
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(rho(p, T)*R*T);
}


inline Foam::scalar Foam::rPolynomialLaw::psi
(
    const scalar p,
    const scalar T
) const
{
    return sqr(rho(p, T))*(C_[3] + C_[4]*T);
}


inline Foam::scalar Foam::rPolynomialLaw::CpMCv
(
    const scalar p,
    const scalar T
) const
{
    return T*sqr(C_[1] + 2*C_[2]*T - C_[4]*p)/(C_[3] + C_[4]*T);
}


inline Foam::scalar Foam::rPolynomialLaw::alphav
(
    const scalar p,
    const scalar T
) const
{
    return rho(p, T)*(C_[1] + 2*C_[2]*T - C_[4]*p);
}

// ************************************************************************* //
