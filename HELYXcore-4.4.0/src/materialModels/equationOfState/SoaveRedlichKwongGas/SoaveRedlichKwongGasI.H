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

#include "SoaveRedlichKwongGas.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::SoaveRedlichKwongGas>
Foam::SoaveRedlichKwongGas::clone() const
{
    return autoPtr<SoaveRedlichKwongGas>(new SoaveRedlichKwongGas(*this));
}


inline Foam::autoPtr<Foam::SoaveRedlichKwongGas>
Foam::SoaveRedlichKwongGas::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<SoaveRedlichKwongGas>
    (
        new SoaveRedlichKwongGas(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::SoaveRedlichKwongGas::alpha(scalar T) const
{
    return sqr(1 + kappa_*(1 - sqrt(T/Tc_)));
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::ap(scalar T) const
{
    return kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T*Tc_));
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::app(scalar T) const
{
    return kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T)*Tc_));
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::A(scalar p, scalar T) const
{
    return a_*alpha(T)*p/sqr(RR*T);
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::B(scalar p, scalar T) const
{
    return b_*p/(RR*T);
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::rho
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::hContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = Z(p, T)*RR*T/p;

    return p*v - RR*T + 1/b_*log(1 + b_/v)*(a_*alpha(T) - T*ap(T));
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::CpContribution
(
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);
    const scalar Z = this->Z(p, T);
    const scalar ap = this->ap(T);
    const scalar app = this->app(T);

    const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return
    (
        app*(T/b_)*log(1 + B/Z) + RR*sqr(M - N)/(sqr(M) - A*(2*Z + B)) - RR
    )/W;
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::eContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = Z(p, T)*RR*T/p;

    return 1/b_*log(1 + b_/v)*(a_*alpha(T) - T*ap(T));
}



inline Foam::scalar Foam::SoaveRedlichKwongGas::CvContribution
(
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    const scalar B = this->B(p,T);
    const scalar Z = this->Z(p, T);
    const scalar app = this->app(T);

    return (app*(T/b_)*log(1 + B/Z) - RR)/W;
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::spContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = Z(p, T)*RR*T/p;

    return RR*log(1 - b_/v) - 1/b_*log(1 + b_/v)*T*ap(T);
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::svContribution() const
{
    NotImplemented;
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::psi
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return 1.0/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::Z
(
    const scalar p,
    const scalar T
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);

    const scalar a2 = - 1;
    const scalar a1 = A - B - B*B;
    const scalar a0 = -B*A;

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = Foam::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*Foam::cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + constant::mathematical::twoPi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = cbrt(Rl + D05);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -cbrt(mag(Rl - D05));
        }
        else
        {
            Tl = cbrt(Rl - D05);
        }

        root = S + Tl - a2/3.0;
    }
    return root;
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::CpMCv
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);
    const scalar Z = this->Z(p, T);
    const scalar ap = this->ap(T);

    const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return R*sqr(M - N)/(sqr(M) - A*(2*Z + B));
}


inline Foam::scalar Foam::SoaveRedlichKwongGas::alphav() const
{
    NotImplemented;
}


// ************************************************************************* //
