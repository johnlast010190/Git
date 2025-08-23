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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/icoPolynomial/icoPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoPolynomial<Specie, PolySize>::icoPolynomial
(
    const Specie& sp,
    const fvMesh& mesh,
    const Polynomial<PolySize>& rhoCoeffs,
    scalar rhoMin,
    scalar rhoMax
)
:
    Specie(sp),
    rhoCoeffs_(rhoCoeffs),
    rhoMin_(rhoMin),
    rhoMax_(rhoMax)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoPolynomial<Specie, PolySize>::icoPolynomial
(
    const word& name,
    const icoPolynomial<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_),
    rhoMin_(ip.rhoMin_),
    rhoMax_(ip.rhoMax_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<Specie, PolySize>>
Foam::icoPolynomial<Specie, PolySize>::clone() const
{
    return autoPtr<icoPolynomial<Specie, PolySize>>
    (
        new icoPolynomial<Specie, PolySize>(*this)
    );
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<Specie, PolySize>>
Foam::icoPolynomial<Specie, PolySize>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<icoPolynomial<Specie, PolySize>>
    (
        new icoPolynomial<Specie, PolySize>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::rho
(
    scalar p,
    scalar T
) const
{
    return max(rhoMin_, min(rhoMax_, rhoCoeffs_.value(T)));
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::h
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::e
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::sp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::sv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::Z
(
    scalar p,
    scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::CpMCv
(
    scalar p,
    scalar T
) const
{
    scalar drdt = rhoCoeffs_.derivative(T);
    scalar rw = rhoCoeffs_.value(T);

    if (rw > rhoMax_)
    {
        rw = rhoMax_;
        drdt = 0;
    }
    else if (rw < rhoMin_)
    {
        rw = rhoMin_;
        drdt = 0;
    }

    return -(p/sqr(rw))*drdt;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::alphav
(
    scalar p,
    scalar T
) const
{
    return -rhoCoeffs_.derivative(T)/rho(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline void Foam::icoPolynomial<Specie, PolySize>::operator+=
(
    const icoPolynomial<Specie, PolySize>& ip
)
{
    NotImplemented;
}


template<class Specie, int PolySize>
inline void Foam::icoPolynomial<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator+
(
    const icoPolynomial<Specie, PolySize>& ip1,
    const icoPolynomial<Specie, PolySize>& ip2
)
{
    NotImplemented;
}


template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const icoPolynomial<Specie, PolySize>& ip
)
{
    return icoPolynomial<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.mesh(),
        ip.rhoCoeffs_,
        ip.rhoMin_,
        ip.rhoMax_
    );
}


template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator==
(
    const icoPolynomial<Specie, PolySize>& ip1,
    const icoPolynomial<Specie, PolySize>& ip2
)
{
    NotImplemented;
}


// ************************************************************************* //
