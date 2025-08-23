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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::ePolynomial> Foam::ePolynomial::clone() const
{
    return autoPtr<ePolynomial>(new ePolynomial(*this));
}


inline Foam::autoPtr<Foam::ePolynomial> Foam::ePolynomial::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<ePolynomial>
    (
        new ePolynomial(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ePolynomial::Cv
(
    const scalar T,
    const scalar CvContribution
) const
{
    return CvContribution + CvCoeffs_.value(T);
}


inline Foam::scalar Foam::ePolynomial::es
(
    const scalar T,
    const scalar eContribution
) const
{
    return eContribution + eCoeffs_.value(T);
}


inline Foam::scalar Foam::ePolynomial::ea
(
    const scalar T,
    const scalar eContribution
) const
{
    return es(T, eContribution) + hf();
}


inline Foam::scalar Foam::ePolynomial::Cp
(
    const scalar Cv,
    const scalar CpMCv
) const
{
    return Cv + CpMCv;
}


inline Foam::scalar Foam::ePolynomial::hs
(
    const scalar p,
    const scalar es,
    const scalar rho
) const
{
    return es + p/rho;
}


inline Foam::scalar Foam::ePolynomial::ha
(
    const scalar p,
    const scalar ea,
    const scalar rho
) const
{
    return ea + p/rho;
}


inline Foam::scalar Foam::ePolynomial::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::ePolynomial::s
(
    const scalar T,
    const scalar svContribution
) const
{
    return svContribution + eCoeffs_.value(T);
}


inline Foam::scalar Foam::ePolynomial::dCpdT(const scalar T) const
{
    NotImplemented;
    return CvCoeffs_.derivative(T); // + EquationOfState::dCpdT
}

// inline Foam::scalar Foam::ePolynomial::gStd(const scalar T) const
// {
//     return hCoeffs_.value(T) - sCoeffs_.value(T)*T;
//  + Pstd/rho(Pstd, T); // This requires coeff mixtures for EOS
// }


// inline Foam::scalar Foam::ePolynomial::deltaS(const scalar T) const
// {
//     //TODO: Introduce coefficient mixture of equation of state
//     return sCoeffs_.value(T);
// }


// ************************************************************************* //
