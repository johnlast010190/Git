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

inline Foam::autoPtr<Foam::eConst> Foam::eConst::clone() const
{
    return autoPtr<eConst>(new eConst(*this));
}


inline Foam::autoPtr<Foam::eConst> Foam::eConst::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<eConst>(new eConst(obr, dict, phaseName, specieName, name));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::eConst::Cv(const scalar CvContribution) const
{
    return CvContribution + Cv_;
}


inline Foam::scalar Foam::eConst::es
(
    const scalar T,
    const scalar eContribution
) const
{
    return Cv_*(T - Tref_) + esref_ + eContribution;
}


inline Foam::scalar Foam::eConst::ea
(
    const scalar T,
    const scalar eContribution
) const
{
    return es(T, eContribution) + hf();
}


inline Foam::scalar Foam::eConst::Cp(const scalar Cv, const scalar CpMCv) const
{
    return Cv + CpMCv;
}


inline Foam::scalar Foam::eConst::hs
(
    const scalar p,
    const scalar es,
    const scalar rho
) const
{
    return es + p/rho;
}


inline Foam::scalar Foam::eConst::ha
(
    const scalar p,
    const scalar ea,
    const scalar rho
) const
{
    return ea + p/rho;
}


inline Foam::scalar Foam::eConst::s
(
    const scalar T,
    const scalar svContribution,
    const scalar Cp
) const
{
    return Cp*log(T/Tstd) + svContribution;
}


inline Foam::scalar Foam::eConst::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::eConst::dCpdT(const scalar T) const
{
    NotImplemented;
}

// We need mixtures of the rho models in order to implement the "e" based
// chemistry coefficient mixtures.
// inline Foam::scalar Foam::eConst::gStd(const scalar T) const
// {
//     return
//         Cv_*(T - Tref_) + esref_ + hf() // + Pstd/EquationOfState::rho(Pstd, T)
//       - S(Pstd, T)*T;
// }



// ************************************************************************* //
