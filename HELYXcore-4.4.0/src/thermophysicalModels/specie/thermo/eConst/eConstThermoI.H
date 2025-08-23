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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const EquationOfState& st,
    const scalar Cv,
    const scalar hf,
    const scalar Tref,
    const scalar esref
)
:
    EquationOfState(st),
    Cv_(Cv),
    hf_(hf),
    Tref_(Tref),
    esref_(esref)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const word& name,
    const eConstThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    hf_(ct.hf_),
    Tref_(ct.Tref_),
    esref_(ct.esref_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState>>
Foam::eConstThermo<EquationOfState>::clone() const
{
    return autoPtr<eConstThermo<EquationOfState>>
    (
        new eConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState>>
Foam::eConstThermo<EquationOfState>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<eConstThermo<EquationOfState>>
    (
        new eConstThermo<EquationOfState>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::es
(
    const scalar p,
    const scalar T
) const
{
    return Cv_*(T - Tref_) + esref_ + EquationOfState::e(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::hf() const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T)*log(T/Tstd) + EquationOfState::sv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return
        Cv_*(T - Tref_) + esref_ + hf() + Pstd/EquationOfState::rho(Pstd, T)
      - s(Pstd, T)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eConstThermo<EquationOfState>::operator+=
(
    const eConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        if
        (
            eConstThermo<EquationOfState>::debug
         && notEqual(Tref_, ct.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << Tref_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << ct.Tref_ << " for "
                << (ct.name().size() ? ct.name() : "others")
                << exit(FatalError);
        }
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        hf_ = Y1*hf_ + Y2*ct.hf_;
        esref_ = Y1*esref_ + Y2*ct.esref_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator+
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return eConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.hf_,
            ct1.Tref_,
            ct1.esref_
        );
    }
    else
    {
        if
        (
            eConstThermo<EquationOfState>::debug
         && notEqual(ct1.Tref_, ct2.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << ct1.Tref_ << " for "
                << (ct1.name().size() ? ct1.name() : "others")
                << " != " << ct2.Tref_ << " for "
                << (ct2.name().size() ? ct2.name() : "others")
                << exit(FatalError);
        }
        return eConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.hf_
          + ct2.Y()/eofs.Y()*ct2.hf_,
            ct1.Tref_,
            ct1.Y()/eofs.Y()*ct1.esref_
          + ct2.Y()/eofs.Y()*ct2.esref_
        );
    }
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eConstThermo<EquationOfState>& ct
)
{
    return eConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.hf_,
        ct.Tref_,
        ct.esref_
    );
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator==
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    if
    (
        eConstThermo<EquationOfState>::debug
     && notEqual(ct1.Tref_, ct2.Tref_)
    )
    {
        FatalErrorInFunction
            << "Tref " << ct1.Tref_ << " for "
            << (ct1.name().size() ? ct1.name() : "others")
            << " != " << ct2.Tref_ << " for "
            << (ct2.name().size() ? ct2.name() : "others")
            << exit(FatalError);
    }

    return eConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.hf_
      - ct1.Y()/eofs.Y()*ct1.hf_,
        ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.esref_
      - ct1.Y()/eofs.Y()*ct1.esref_
    );
}


// ************************************************************************* //
