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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermo/hPolynomial/hPolynomialThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const EquationOfState& pt,
    const scalar hf,
    const scalar Sf,
    const Polynomial<PolySize>& CpCoeffs,
    const typename Polynomial<PolySize>::intPolyType& hCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    EquationOfState(pt),
    hf_(hf),
    Sf_(Sf),
    CpCoeffs_(CpCoeffs),
    hCoeffs_(hCoeffs),
    sCoeffs_(sCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const word& name,
    const hPolynomialThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    Sf_(pt.Sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::hs
(
    const scalar p,
    const scalar T
) const
{
    return hCoeffs_.value(T) + EquationOfState::h(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::hf()
const
{
    return hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::s
(
    const scalar p,
    const scalar T
) const
{
    return sCoeffs_.value(T) + EquationOfState::sp(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::gStd
(
    const scalar T
) const
{
    return hCoeffs_.value(T) - sCoeffs_.value(T)*T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.derivative(T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator=
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    EquationOfState::operator=(pt);

    hf_ = pt.hf_;
    Sf_ = pt.Sf_;
    CpCoeffs_ = pt.CpCoeffs_;
    hCoeffs_ = pt.hCoeffs_;
    sCoeffs_ = pt.sCoeffs_;
}


template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator+=
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pt.Y()/this->Y();

        hf_ = Y1*hf_ + Y2*pt.hf_;
        Sf_ = Y1*Sf_ + Y2*pt.Sf_;
        CpCoeffs_ = Y1*CpCoeffs_ + Y2*pt.CpCoeffs_;
        hCoeffs_ = Y1*hCoeffs_ + Y2*pt.hCoeffs_;
        sCoeffs_ = Y1*sCoeffs_ + Y2*pt.sCoeffs_;
    }
}


template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            pt1.hf_,
            pt1.Sf_,
            pt1.CpCoeffs_,
            pt1.hCoeffs_,
            pt1.sCoeffs_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            Y1*pt1.hf_ + Y2*pt2.hf_,
            Y1*pt1.Sf_ + Y2*pt2.Sf_,
            Y1*pt1.CpCoeffs_ + Y2*pt2.CpCoeffs_,
            Y1*pt1.hCoeffs_ + Y2*pt2.hCoeffs_,
            Y1*pt1.sCoeffs_ + Y2*pt2.sCoeffs_
        );
    }
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    return hPolynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.hf_,
        pt.Sf_,
        pt.CpCoeffs_,
        pt.hCoeffs_,
        pt.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return hPolynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        Y2*pt2.hf_       - Y1*pt1.hf_,
        Y2*pt2.Sf_       - Y1*pt1.Sf_,
        Y2*pt2.CpCoeffs_ - Y1*pt1.CpCoeffs_,
        Y2*pt2.hCoeffs_  - Y1*pt1.hCoeffs_,
        Y2*pt2.sCoeffs_  - Y1*pt1.sCoeffs_
    );
}


// ************************************************************************* //
