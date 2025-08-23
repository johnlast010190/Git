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
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "thermo/hTabulated/hTabulatedThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const EquationOfState& st,
    const dictionary& CpDict,
    const autoPtr<interpolation2DTable<scalar>>& CpTable,
    const scalar hf
)
:
    EquationOfState(st),
    CpDict_(CpDict),
    CpTable_(new interpolation2DTable<scalar>(CpTable())),
    hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const word& name,
    const hTabulatedThermo& ct
)
:
    EquationOfState(name, ct),
    CpDict_(ct.CpDict_),
    CpTable_(new interpolation2DTable<scalar>(ct.CpTable_())),
    hf_(ct.hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hTabulatedThermo<EquationOfState>>
Foam::hTabulatedThermo<EquationOfState>::clone() const
{
    return autoPtr< hTabulatedThermo<EquationOfState>>
    (
        new  hTabulatedThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hTabulatedThermo<EquationOfState>>
Foam::hTabulatedThermo<EquationOfState>::New(const objectRegistry& obr, const dictionary& dict)
{
    return autoPtr< hTabulatedThermo<EquationOfState>>
    (
        new  hTabulatedThermo<EquationOfState>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_()(p,T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_()(p,T)*T + hf_ + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_()(p,T)*T + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::hf() const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::s
(
    const scalar p, const scalar T
) const
{
    return CpTable_()(p,T)*log(T/Tstd) + EquationOfState::sp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hTabulatedThermo<EquationOfState>::operator+=
(
    const  hTabulatedThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        CpTable_().idaxpy(Y1, Y2, ct.CpTable_());
        hf_ = Y1*hf_ + Y2*ct.hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator+
(
    const  hTabulatedThermo<EquationOfState>& ct1,
    const  hTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return hTabulatedThermo<EquationOfState>
        (
            eofs,
            ct1.CpDict_,
            ct1.CpTable_,
            ct1.hf_
        );
    }
    else
    {
        interpolation2DTable<scalar> iCpT(ct1.CpTable_());
        iCpT.idaxpy
        (
            ct1.Y()/eofs.Y(), ct2.Y()/eofs.Y(), ct2.CpTable_()
        );


        return hTabulatedThermo<EquationOfState>
        (
            eofs,
            dictionary("emptyCpDictPlus"),
            iCpT,
            ct1.Y()/eofs.Y()*ct1.hf_
          + ct2.Y()/eofs.Y()*ct2.hf_
        );
    }
}


template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const  hTabulatedThermo<EquationOfState>& ct
)
{
    return  hTabulatedThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.CpDict_,
        ct.CpTable_,
        ct.hf_
    );
}


template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator==
(
    const  hTabulatedThermo<EquationOfState>& ct1,
    const  hTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    autoPtr<interpolation2DTable<scalar>> iCpT
    (
        new interpolation2DTable<scalar>(ct1.CpTable_())
    );
    iCpT->idaxpy
    (
        ct1.Y()/eofs.Y(),
        -ct2.Y()/eofs.Y(),
        ct2.CpTable_()
    );

    return  hTabulatedThermo<EquationOfState>
    (
        eofs,
        dictionary("emptyCpDictMinus"),
        iCpT,
        ct1.Y()/eofs.Y()*ct1.hf_
      - ct2.Y()/eofs.Y()*ct2.hf_
    );
}


// ************************************************************************* //
