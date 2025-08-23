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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hFunctionT> Foam::hFunctionT::clone() const
{
    return autoPtr<hFunctionT>(new hFunctionT(*this));
}


inline Foam::autoPtr<Foam::hFunctionT> Foam::hFunctionT::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hFunctionT>
    (
        new hFunctionT(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hFunctionT::Cp
(
    const scalar T,
    const scalar CpContribution
) const
{
    return CpContribution + Cp_().value(T);
}


inline Foam::scalar Foam::hFunctionT::ha
(
    const scalar T,
    const scalar hContribution
) const
{
    return
        hContribution
      + (h_.valid() ? h_->value(T) : (Cp_().integrate(0, T) + hf_));
}


inline Foam::scalar Foam::hFunctionT::hs
(
    const scalar ha,
    const scalar hf
) const
{
    return ha - hf;
}


inline Foam::scalar Foam::hFunctionT::s
(
    const scalar T,
    const scalar spContribution
) const
{
    return spContribution + Cp_().integrateYoverX(Tstd, T) + Sf_;
}


inline Foam::scalar Foam::hFunctionT::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::hFunctionT::dCpdT(const scalar T) const
{
    NotImplemented;
}


inline Foam::scalar Foam::hFunctionT::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::hFunctionT::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::hFunctionT::ea
(
    const scalar p,
    const scalar rho,
    const scalar Ha
) const
{
    return Ha - p/rho;
}


// ************************************************************************* //
