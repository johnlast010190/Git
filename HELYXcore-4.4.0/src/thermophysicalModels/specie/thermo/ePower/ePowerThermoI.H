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

#include "thermo/ePower/ePowerThermo.H"
#include "specie/specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::ePowerThermo<EquationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0)
    {
        FatalErrorInFunction
            << "attempt to evaluate ePowerThermo<EquationOfState>"
               " for negative temperature " << T
            << abort(FatalError);
    }
}


template<class EquationOfState>
inline Foam::ePowerThermo<EquationOfState>::ePowerThermo
(
    const word& name,
    const ePowerThermo& jt
)
:
    EquationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    hf_(jt.hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::ePowerThermo<EquationOfState>::ePowerThermo
(
    const EquationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar hf
)
:
    EquationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    hf_(hf)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::ePowerThermo<EquationOfState>>
Foam::ePowerThermo<EquationOfState>::clone() const
{
    return autoPtr<ePowerThermo<EquationOfState>>
    (
        new ePowerThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::ePowerThermo<EquationOfState>>
Foam::ePowerThermo<EquationOfState>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<ePowerThermo<EquationOfState>>
    (
        new ePowerThermo<EquationOfState>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_) + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::es
(
    const scalar p,
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + EquationOfState::e(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::hf() const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_)
      + EquationOfState::sv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + hf() + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::ePowerThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::ePowerThermo<EquationOfState>::operator+=
(
    const ePowerThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        hf_ = Y1*hf_ + Y2*ct.hf_;
        c0_ = Y1*c0_ + Y2*ct.c0_;
        n0_ = Y1*n0_ + Y2*ct.n0_;
        Tref_ = Y1*Tref_ + Y2*ct.Tref_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::ePowerThermo<EquationOfState> Foam::operator+
(
    const ePowerThermo<EquationOfState>& ct1,
    const ePowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return ePowerThermo<EquationOfState>
        (
            eofs,
            ct1.c0_,
            ct1.n0_,
            ct1.Tref_,
            ct1.hf_
        );
    }
    else
    {
        return ePowerThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.c0_
          + ct2.Y()/eofs.Y()*ct2.c0_,
            ct1.Y()/eofs.Y()*ct1.n0_
          + ct2.Y()/eofs.Y()*ct2.n0_,
            ct1.Y()/eofs.Y()*ct1.Tref_
          + ct2.Y()/eofs.Y()*ct2.Tref_,
            ct1.Y()/eofs.Y()*ct1.hf_
          + ct2.Y()/eofs.Y()*ct2.hf_
        );
    }
}


template<class EquationOfState>
inline Foam::ePowerThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const ePowerThermo<EquationOfState>& ct
)
{
    return ePowerThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.hf_
    );
}


template<class EquationOfState>
inline Foam::ePowerThermo<EquationOfState> Foam::operator==
(
    const ePowerThermo<EquationOfState>& ct1,
    const ePowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return ePowerThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.c0_
      - ct1.Y()/eofs.Y()*ct1.c0_,
        ct2.Y()/eofs.Y()*ct2.n0_
      - ct1.Y()/eofs.Y()*ct1.n0_,
        ct2.Y()/eofs.Y()*ct2.Tref_
      - ct1.Y()/eofs.Y()*ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.hf_
      - ct1.Y()/eofs.Y()*ct1.hf_
    );
}


// ************************************************************************* //