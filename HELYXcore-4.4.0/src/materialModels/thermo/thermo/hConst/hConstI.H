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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hConst> Foam::hConst::clone() const
{
    return autoPtr<hConst>(new hConst(*this));
}


inline Foam::autoPtr<Foam::hConst> Foam::hConst::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hConst>(new hConst(obr, dict, phaseName, specieName, name));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hConst::Cp(const scalar CPContribution) const
{
    return CPContribution + Cp_;
}


inline Foam::scalar Foam::hConst::ha
(
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + Cp_*(T - Tref_) + hsref_ + hf_;
}


inline Foam::scalar Foam::hConst::hs
(
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + Cp_*(T - Tref_) + hsref_;
}


inline Foam::scalar Foam::hConst::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::hConst::s
(
    const scalar T,
    const scalar sContribution
) const
{
    return sContribution + Sf_ + Cp_*log(T/Tstd);
}


inline Foam::scalar Foam::hConst::dCpdT(const scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::hConst::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::hConst::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::hConst::ea
(
    const scalar p,
    const scalar rho,
    const scalar Ha
) const
{
    return Ha - p/rho;
}


inline Foam::scalar Foam::hConst::gStd(const scalar T) const
{
    return Cp_*(T - Tref_) + hsref_ + hf_ - Cp_*T*log(T/Tstd);
}


inline Foam::scalar Foam::hConst::deltaS(const scalar T) const
{
    //TODO: Introduce coefficient mixture of equation of state
    return Sf_ + Cp_*log(T/Tstd);
}


// ************************************************************************* //
