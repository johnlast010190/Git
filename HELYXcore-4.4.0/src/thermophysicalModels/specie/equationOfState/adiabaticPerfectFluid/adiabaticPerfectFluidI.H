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
    (c) 2016 Engys Ltd.
    (c) 2013-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const Specie& sp,
    const fvMesh& mesh,
    const scalar p0,
    const scalar rho0,
    const scalar gamma,
    const scalar B
)
:
    Specie(sp),
    p0_(p0),
    rho0_(rho0),
    gamma_(gamma),
    B_(B)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const word& name,
    const adiabaticPerfectFluid<Specie>& pf
)
:
    Specie(name, pf),
    p0_(pf.p0_),
    rho0_(pf.rho0_),
    gamma_(pf.gamma_),
    B_(pf.B_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::adiabaticPerfectFluid<Specie>>
Foam::adiabaticPerfectFluid<Specie>::clone() const
{
    return autoPtr<adiabaticPerfectFluid<Specie>>
    (
        new adiabaticPerfectFluid<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::adiabaticPerfectFluid<Specie>>
Foam::adiabaticPerfectFluid<Specie>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<adiabaticPerfectFluid<Specie>>
    (
        new adiabaticPerfectFluid<Specie>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho0_*pow((p + B_)/(p0_ + B_), 1.0/gamma_);
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::h
(
    scalar p,
    scalar T
) const
{
    // ***HGW This needs to be added, it is not 0
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::e
(
    scalar p,
    scalar T
) const
{
    // ***HGW This needs to be added, it is H - p/rho
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::sp
(
    scalar p,
    scalar T
) const
{
    scalar n = 1 - 1.0/gamma_;
    return
       -pow(p0_ + B_, 1.0/gamma_)*(pow((p + B_), n) - pow((Pstd + B_), n))
       /(rho0_*T*n);
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::sv
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return
        (rho0_/(gamma_*(p0_ + B_)))*pow((p + B_)/(p0_ + B_), 1.0/gamma_ - 1.0);
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::adiabaticPerfectFluid<Specie>::alphav
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::adiabaticPerfectFluid<Specie>::operator+=
(
    const adiabaticPerfectFluid<Specie>& pf
)
{
    NotImplemented;
}


template<class Specie>
inline void Foam::adiabaticPerfectFluid<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::adiabaticPerfectFluid<Specie> Foam::operator+
(
    const adiabaticPerfectFluid<Specie>& pf1,
    const adiabaticPerfectFluid<Specie>& pf2
)
{
    NotImplemented;
}


template<class Specie>
inline Foam::adiabaticPerfectFluid<Specie> Foam::operator*
(
    const scalar s,
    const adiabaticPerfectFluid<Specie>& pf
)
{
    return adiabaticPerfectFluid<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.mesh(),
        pf.p0_,
        pf.rho0_,
        pf.gamma_,
        pf.B_
    );
}


template<class Specie>
inline Foam::adiabaticPerfectFluid<Specie> Foam::operator==
(
    const adiabaticPerfectFluid<Specie>& pf1,
    const adiabaticPerfectFluid<Specie>& pf2
)
{
    NotImplemented;
}


// ************************************************************************* //
