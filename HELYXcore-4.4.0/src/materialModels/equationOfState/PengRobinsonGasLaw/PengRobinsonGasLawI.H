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
    (c) 2014-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "PengRobinsonGasLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::PengRobinsonGasLaw>
Foam::PengRobinsonGasLaw::clone() const
{
    return autoPtr<PengRobinsonGasLaw>(new PengRobinsonGasLaw(*this));
}


inline Foam::autoPtr<Foam::PengRobinsonGasLaw> Foam::PengRobinsonGasLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<PengRobinsonGasLaw>
    (
        new PengRobinsonGasLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::PengRobinsonGasLaw::alpha(scalar T) const
{
    return sqr(1 + kappa_*(1 - sqrt(T/Tc_)));
}


inline Foam::scalar Foam::PengRobinsonGasLaw::B(scalar p, scalar T) const
{
    return b_*p/(RR*T);
}


inline Foam::scalar Foam::PengRobinsonGasLaw::A(scalar p, scalar T) const
{
    return a_*alpha(T)*p/sqr(RR*T);
}


inline Foam::scalar Foam::PengRobinsonGasLaw::rho
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return p/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::PengRobinsonGasLaw::hContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar Z = this->Z(p, T);

    return
        R*Tc_
       *(
           Tr*(Z - 1)
         - 2.078*(1 + kappa_)*sqrt(alpha(T))*log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


inline Foam::scalar Foam::PengRobinsonGasLaw::CpContribution
(
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);
    const scalar ap = kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T*Tc_));
    const scalar app = kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T)*Tc_));
    const scalar Z = this->Z(p, T);
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    const scalar root2 = sqrt(2.0);

    return
    (
        app*(T/(2*root2*b_))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
      - RR
    )/W;
}


inline Foam::scalar Foam::PengRobinsonGasLaw::eContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar Pr = p/Pc_;
    const scalar B = 0.07780*Pr/(T/Tc_);
    const scalar Z = this->Z(p, T);

    return
        R*Tc_
       *(
            -2.078*(1 + kappa_)*sqrt(alpha(T))*log((Z + 2.414*B)/(Z - 0.414*B))
        );
}



inline Foam::scalar Foam::PengRobinsonGasLaw::CvContribution
(
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    const scalar B = this->B(p, T);
    const scalar Z = this->Z(p, T);
    const scalar app = kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T)*Tc_));
    const scalar root2 = sqrt(2.0);

    return
    (
        app*(T/(2*root2*b_))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      - RR
    )/W;
}


inline Foam::scalar Foam::PengRobinsonGasLaw::spContribution
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar Z = this->Z(p, T);

    return
        R
       *(
          - log(p/Pstd)
          + (
                log(Z - B)
              - 2.078*kappa_*((1 + kappa_)/sqrt(Tr) - kappa_)
               *log((Z + 2.414*B)/(Z - 0.414*B))
            )
        );
}


inline Foam::scalar Foam::PengRobinsonGasLaw::svContribution() const
{
    NotImplemented;
}


inline Foam::scalar Foam::PengRobinsonGasLaw::psi
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    return 1.0/(Z(p, T)*R*T);
}


inline Foam::scalar Foam::PengRobinsonGasLaw::Z
(
    const scalar p,
    const scalar T
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);

    const scalar a2 = B - 1;
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = Foam::acos(Rl/Foam::sqrt(-Q3));
        const scalar qm = 2*Foam::sqrt(-Q);
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


inline Foam::scalar Foam::PengRobinsonGasLaw::CpMCv
(
    const scalar p,
    const scalar T,
    const scalar R
) const
{
    const scalar A = this->A(p, T);
    const scalar B = this->B(p, T);
    const scalar ap = kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T*Tc_));
    const scalar Z = this->Z(p, T);
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return R*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}


inline Foam::scalar Foam::PengRobinsonGasLaw::alphav() const
{
    NotImplemented;
}


// ************************************************************************* //
