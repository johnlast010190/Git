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
    (c) 2020-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "icoTabulated.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::icoTabulated<Specie>::icoTabulated
(
    const Specie& sp,
    const nonUniformTable& rho
)
:
    Specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::icoTabulated<Specie>::icoTabulated
(
    const word& name,
    const icoTabulated<Specie>& ip
)
:
    Specie(name, ip),
    rho_(ip.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::icoTabulated<Specie>>
Foam::icoTabulated<Specie>::clone() const
{
    return autoPtr<icoTabulated<Specie>>
    (
        new icoTabulated<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::icoTabulated<Specie>>
Foam::icoTabulated<Specie>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<icoTabulated<Specie>>(new icoTabulated<Specie>(obr, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho_.value(T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::h
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::e
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::sp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::sv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::alphav
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
}


// ************************************************************************* //