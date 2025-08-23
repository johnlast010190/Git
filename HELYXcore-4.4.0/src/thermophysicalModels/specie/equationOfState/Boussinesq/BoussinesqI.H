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
    (c) 2015-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/Boussinesq/Boussinesq.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Boussinesq<Specie>::Boussinesq
(
    const Specie& sp,
    const fvMesh& mesh,
    const scalar rho0,
    const scalar T0,
    const scalar beta
)
:
    Specie(sp),
    rho0_(rho0),
    T0_(T0),
    beta_(beta)
{}


template<class Specie>
inline Foam::Boussinesq<Specie>::Boussinesq
(
    const word& name,
    const Boussinesq<Specie>& b
)
:
    Specie(name, b),
    rho0_(b.rho0_),
    T0_(b.T0_),
    beta_(b.beta_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Boussinesq<Specie>>
Foam::Boussinesq<Specie>::clone() const
{
    return autoPtr<Boussinesq<Specie>>(new Boussinesq<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::Boussinesq<Specie>>
Foam::Boussinesq<Specie>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<Boussinesq<Specie>>(new Boussinesq<Specie>(obr, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::rho(scalar p, scalar T) const
{
    return rho0_*(1.0 - beta_*(T - T0_));
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::h(scalar p, scalar T) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::e(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::sp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::sv(scalar p, scalar T) const
{
    NotImplemented;
}



template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::psi(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::alphav(scalar p, scalar T) const
{
    return rho0_/this->rho(p, T)*beta_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Boussinesq<Specie>::operator+=
(
    const Boussinesq<Specie>& b
)
{
    NotImplemented;
}


template<class Specie>
inline void Foam::Boussinesq<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator+
(
    const Boussinesq<Specie>& b1,
    const Boussinesq<Specie>& b2
)
{
    NotImplemented;
}


template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator*
(
    const scalar s,
    const Boussinesq<Specie>& b
)
{
    return Boussinesq<Specie>
    (
        s*static_cast<const Specie&>(b),
        b.mesh(),
        b.rho0_,
        b.T0_,
        b.beta_
    );
}


template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator==
(
    const Boussinesq<Specie>& b1,
    const Boussinesq<Specie>& b2
)
{
    NotImplemented;
}


// ************************************************************************* //
