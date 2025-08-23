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

inline Foam::autoPtr<Foam::hPolynomial> Foam::hPolynomial::clone() const
{
    return autoPtr<hPolynomial>(new hPolynomial(*this));
}


inline Foam::autoPtr<Foam::hPolynomial> Foam::hPolynomial::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hPolynomial>
    (
        new hPolynomial(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hPolynomial::Cp
(
    const scalar T,
    const scalar CpContribution
) const
{
    return CpContribution + CpCoeffs_.value(T);
}


inline Foam::scalar Foam::hPolynomial::ha
(
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + hCoeffs_.value(T) + hf_;
}


inline Foam::scalar Foam::hPolynomial::hs
(
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + hCoeffs_.value(T);
}


inline Foam::scalar Foam::hPolynomial::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::hPolynomial::s
(
    const scalar T,
    const scalar SContribution
) const
{
    return SContribution + sCoeffs_.value(T);
}


inline Foam::scalar Foam::hPolynomial::dCpdT(const scalar T) const
{
    return CpCoeffs_.derivative(T);
}


inline Foam::scalar Foam::hPolynomial::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::hPolynomial::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::hPolynomial::ea
(
    const scalar p,
    const scalar rho,
    const scalar Ha
) const
{
    return Ha - p/rho;
}


inline Foam::scalar Foam::hPolynomial::gStd(const scalar T) const
{
    return hCoeffs_.value(T) - sCoeffs_.value(T)*T;
}


inline Foam::scalar Foam::hPolynomial::deltaS(const scalar T) const
{
    //TODO: Introduce coefficient mixture of equation of state
    return sCoeffs_.value(T);
}


// ************************************************************************* //
