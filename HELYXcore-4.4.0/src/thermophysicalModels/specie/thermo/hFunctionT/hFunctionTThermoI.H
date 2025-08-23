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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "thermo/hFunctionT/hFunctionTThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState>::hFunctionTThermo
(
    const EquationOfState& pt,
    const scalar hf,
    const scalar Sf,
    const Function1<scalar>& Cp
)
:
    EquationOfState(pt),
    hf_(hf),
    Sf_(Sf),
    Cp_(Cp.clone().ptr())
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState>::hFunctionTThermo
(
    const word& name,
    const hFunctionTThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    Sf_(pt.Sf_),
    Cp_(pt.Cp_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{
    return Cp_().value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return Cp_().integrate(0, T) + hf_ + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return Cp_().integrateYoverX(Tstd, T) + Sf_ + EquationOfState::sp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator=
(
    const hFunctionTThermo<EquationOfState>& pt
)
{
    EquationOfState::operator=(pt);

    hf_ = pt.hf_;
    Sf_ = pt.Sf_;
    Cp_.reset(pt.Cp_.clone());
}


template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator+=
(
    const hFunctionTThermo<EquationOfState>& pt
)
{
    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        NotImplemented;
    }
}


template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator+
(
    const hFunctionTThermo<EquationOfState>& pt1,
    const hFunctionTThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return hFunctionTThermo<EquationOfState>
        (
            eofs,
            pt1.hf_,
            pt1.Sf_,
            pt1.Cp_(),
            pt1.h_(),
            pt1.s_()
        );
    }
    {
        NotImplemented;
    }
}


template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hFunctionTThermo<EquationOfState>& pt
)
{
    return hFunctionTThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.hf_,
        pt.Sf_,
        pt.Cp_(),
        pt.h_(),
        pt.s_()
    );
}


template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator==
(
    const hFunctionTThermo<EquationOfState>& pt1,
    const hFunctionTThermo<EquationOfState>& pt2
)
{
    NotImplemented;
}


// ************************************************************************* //
