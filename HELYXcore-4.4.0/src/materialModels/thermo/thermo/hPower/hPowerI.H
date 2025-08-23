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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hPower> Foam::hPower::clone() const
{
    return autoPtr<hPower>(new hPower(*this));
}


inline Foam::autoPtr<Foam::hPower> Foam::hPower::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hPower>(new hPower(obr, dict, phaseName, specieName, name));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hPower::Cp
(
    const scalar T,
    const scalar CpContribution
) const
{
    return CpContribution + c0_*pow(T/Tref_, n0_);
}


inline Foam::scalar Foam::hPower::ha
(
    const scalar hs,
    const scalar hf
) const
{
    return hs + hf;
}


inline Foam::scalar Foam::hPower::hs
(
    const scalar T,
    const scalar hContribution
) const
{
    return
        hContribution
      + c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1));
}


inline Foam::scalar Foam::hPower::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::hPower::s
(
    const scalar T,
    const scalar SContribution
) const
{
    return
        SContribution
      + c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_);
}


inline Foam::scalar Foam::hPower::dCpdT(const scalar T) const
{
    NotImplemented;
}


inline Foam::scalar Foam::hPower::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::hPower::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::hPower::ea
(
    const scalar p,
    const scalar rho,
    const scalar ha
) const
{
    return ha - p/rho;
}


inline Foam::scalar Foam::hPower::gStd(const scalar T) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + hf_
      - c0_*(pow(T, n0_) - pow(Tstd, n0_))*T/(pow(Tref_, n0_)*n0_);
}


inline Foam::scalar Foam::hPower::deltaS(const scalar T) const
{
    NotImplemented;
}


// ************************************************************************* //
