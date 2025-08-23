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

#include "RedlichKwongPengRobinsonGas.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::RedlichKwongPengRobinsonGas>
Foam::RedlichKwongPengRobinsonGas::clone() const
{
    return autoPtr<RedlichKwongPengRobinsonGas>
    (
        new RedlichKwongPengRobinsonGas(*this)
    );
}


inline Foam::autoPtr<Foam::RedlichKwongPengRobinsonGas>
Foam::RedlichKwongPengRobinsonGas::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<RedlichKwongPengRobinsonGas>
    (
        new RedlichKwongPengRobinsonGas(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::alpha(scalar T) const
{
    return pow(3.0/(2 + T/Tc_), kappa_);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::ap(scalar T) const
{
    return -a_*pow(3.0, kappa_)*kappa_/(Tc_*pow(2 + T/Tc_, kappa_ + 1));
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::app(scalar T) const
{
    return
        a_
       *(
            pow(3.0, kappa_)*kappa_*(kappa_ + 1)
           /(sqr(Tc_)*pow(2 + T/Tc_, kappa_ + 2))
        );
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::A
(
    scalar p,
    scalar T
) const
{
    return a_*alpha(T)*p/sqr(RR*T);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::B
(
    scalar p,
    scalar T
) const
{
    return b_*p/(RR*T);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::rho
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::hContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = Z(p, T)*RR*T/p;

    return
        p*v - RR*T
      + 1/(b_*d1MinusD2_)*log((v + d1_*b_)/(v + d2_*b_))
       *(a_*alpha(T) - T*ap(T));
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::CpContribution
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

    const scalar M = (sqr(Z) + d1PlusD2_*B*Z + d1TimesD2_*sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return
    (
        app*(T/(d1MinusD2_*b_))*log((Z + d1_*B)/(Z + d2_*B))
      + RR*sqr(M - N)/(sqr(M) - A*(2*Z + d1PlusD2_*B))
      - RR
    )/W;
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::eContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = Z(p, T)*RR*T/p;

    return
        1/(b_*d1MinusD2_)*log((v + d1_*b_)/(v + d2_*b_))
       *(a_*alpha(T) - T*ap(T));
}



inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::CvContribution
(
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    const scalar B = this->B(p,T);
    const scalar Z = this->Z(p, T);
    const scalar app = this->app(T);

    return
    (
        app*(T/(d1MinusD2_*b_))*log((Z + d1_*B)/(Z  + d2_*B))
      - RR
    )/W;
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::spContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar v = this->Z(p, T)*RR*T/p;

    return
        RR*log(1 - b_/v)
      - 1/(b_*d1MinusD2_)*log((v + d1_*b_)/(v + d2_*b_))*T*ap(T);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::svContribution() const
{
    NotImplemented;
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::psi
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return 1.0/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::Z
(
    const scalar p,
    const scalar T
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);

    const scalar a2 = B*(d1PlusD2_ - 1) - 1;
    const scalar a1 = A + B*(d1TimesD2_*B - d1PlusD2_*(B + 1));
    const scalar a0 = -B*(d1TimesD2_*(sqr(B) + B) + A);

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
            qm*Foam::cos((th + constant::mathematical::twoPi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*Foam::cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

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


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::CpMCv
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

    const scalar M = (sqr(Z) + d1PlusD2_*B*Z + d1TimesD2_*sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return R*sqr(M - N)/(sqr(M) - A*(2*Z + d1PlusD2_*B));
}


inline Foam::scalar Foam::RedlichKwongPengRobinsonGas::alphav() const
{
    NotImplemented;
}


// ************************************************************************* //
