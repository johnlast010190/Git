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

inline Foam::autoPtr<Foam::ePower> Foam::ePower::clone() const
{
    return autoPtr<ePower>(new ePower(*this));
}


inline Foam::autoPtr<Foam::ePower> Foam::ePower::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<ePower>(new ePower(obr, dict, phaseName, specieName, name));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline Foam::scalar Foam::ePower::Cv
(
    const scalar T,
    const scalar CvContribution
) const
{
    return c0_*pow(T/Tref_, n0_) + CvContribution;
}


inline Foam::scalar Foam::ePower::es
(
    const scalar T,
    const scalar eContribution
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + eContribution;
}


inline Foam::scalar Foam::ePower::ea
(
    const scalar T,
    const scalar eContribution
) const
{
    return es(T, eContribution) + hf();
}


inline Foam::scalar Foam::ePower::Cp(const scalar Cv, const scalar CpMCv) const
{
    return Cv + CpMCv;
}


inline Foam::scalar Foam::ePower::hs
(
    const scalar p,
    const scalar es,
    const scalar rho
) const
{
    return es + p/rho;
}


inline Foam::scalar Foam::ePower::ha
(
    const scalar p,
    const scalar ea,
    const scalar rho
) const
{
    return ea + p/rho;
}


inline Foam::scalar Foam::ePower::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::ePower::s
(
    const scalar T,
    const scalar svContribution
) const
{
    return
        svContribution + c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_);
}


inline Foam::scalar Foam::ePower::gStd(const scalar T) const
{
    // Needs adding Sv deparature and Pstd/rho
    NotImplemented;
//     return
//         c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
//       + hf_
//       - c0_*(pow(T, n0_) - pow(Tstd, n0_))*T/(pow(Tref_, n0_)*n0_);
}


inline Foam::scalar Foam::ePower::dCpdT(const scalar T) const
{
    NotImplemented;
}


inline Foam::scalar Foam::ePower::deltaS(const scalar T) const
{
    NotImplemented;
}


// ************************************************************************* //
