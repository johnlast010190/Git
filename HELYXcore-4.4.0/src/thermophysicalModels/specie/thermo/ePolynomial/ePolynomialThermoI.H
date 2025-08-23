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

#include "thermo/ePolynomial/ePolynomialThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const EquationOfState& pt,
    const scalar hf,
    const scalar Sf,
    const Polynomial<PolySize>& CvCoeffs,
    const typename Polynomial<PolySize>::intPolyType& eCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    EquationOfState(pt),
    hf_(hf),
    Sf_(Sf),
    CvCoeffs_(CvCoeffs),
    eCoeffs_(eCoeffs),
    sCoeffs_(sCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const word& name,
    const ePolynomialThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    Sf_(pt.Sf_),
    CvCoeffs_(pt.CvCoeffs_),
    eCoeffs_(pt.eCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return CvCoeffs_.value(T) + EquationOfState::Cv(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::es
(
    const scalar p,
    const scalar T
) const
{
    return eCoeffs_.value(T) + EquationOfState::e(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::hf()
const
{
    return hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::s
(
    const scalar p,
    const scalar T
) const
{
    return sCoeffs_.value(T) + EquationOfState::sv(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::gStd
(
    const scalar T
) const
{
    return
        eCoeffs_.value(T) + hf() + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return CvCoeffs_.derivative(T); // + EquationOfState::dCpdT
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::ePolynomialThermo<EquationOfState, PolySize>::operator+=
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt
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
        CvCoeffs_ = Y1*CvCoeffs_ + Y2*pt.CvCoeffs_;
        eCoeffs_ = Y1*eCoeffs_ + Y2*pt.eCoeffs_;
        sCoeffs_ = Y1*sCoeffs_ + Y2*pt.sCoeffs_;
    }
}


template<class EquationOfState, int PolySize>
inline void Foam::ePolynomialThermo<EquationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt1,
    const ePolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return ePolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            pt1.hf_,
            pt1.Sf_,
            pt1.CvCoeffs_,
            pt1.eCoeffs_,
            pt1.sCoeffs_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return ePolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            Y1*pt1.hf_ + Y2*pt2.hf_,
            Y1*pt1.Sf_ + Y2*pt2.Sf_,
            Y1*pt1.CvCoeffs_ + Y2*pt2.CvCoeffs_,
            Y1*pt1.eCoeffs_ + Y2*pt2.eCoeffs_,
            Y1*pt1.sCoeffs_ + Y2*pt2.sCoeffs_
        );
    }
}


template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const ePolynomialThermo<EquationOfState, PolySize>& pt
)
{
    return ePolynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.hf_,
        pt.Sf_,
        pt.CvCoeffs_,
        pt.eCoeffs_,
        pt.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt1,
    const ePolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return ePolynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        Y2*pt2.hf_       - Y1*pt1.hf_,
        Y2*pt2.Sf_       - Y1*pt1.Sf_,
        Y2*pt2.CvCoeffs_ - Y1*pt1.CvCoeffs_,
        Y2*pt2.eCoeffs_  - Y1*pt1.eCoeffs_,
        Y2*pt2.sCoeffs_  - Y1*pt1.sCoeffs_
    );
}


// ************************************************************************* //