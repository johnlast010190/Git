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

inline Foam::autoPtr<Foam::janaf> Foam::janaf::clone() const
{
    return autoPtr<janaf>(new janaf(*this));
}


inline Foam::autoPtr<Foam::janaf> Foam::janaf::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<janaf>(new janaf(obr, dict, phaseName, specieName, name));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const typename Foam::janaf::coeffArray& Foam::janaf::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


inline Foam::scalar Foam::janaf::limit(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


inline Foam::scalar Foam::janaf::Cp
(
    const scalar T,
    const scalar CpContribution
) const
{
    const coeffArray& a = coeffs(T);
    return CpContribution + ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
}


inline Foam::scalar Foam::janaf::ha
(
    const scalar T,
    const scalar hContribution
) const
{
    const coeffArray& a = coeffs(T);
    return
        hContribution
      + (
            ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
          + a[5]
        );
}


inline Foam::scalar Foam::janaf::hs(const scalar ha, const scalar hf) const
{
    return ha - hf;
}


inline Foam::scalar Foam::janaf::hf() const
{
    const coeffArray& a = lowCpCoeffs_;
    return
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


inline Foam::scalar Foam::janaf::s
(
    const scalar T,
    const scalar SContribution
) const
{
    const coeffArray& a = coeffs(T);
    return
        SContribution
      + (
            (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
          + a[6]
        );
}


inline Foam::scalar Foam::janaf::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::janaf::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::janaf::ea
(
    const scalar p,
    const scalar rho,
    const scalar ha
) const
{
    return ha - p/rho;
}


inline Foam::scalar Foam::janaf::gStd(const scalar T) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (
            a[0]*(1 - log(T))
          - (((a[4]/20.0*T + a[3]/12.0)*T + a[2]/6.0)*T + a[1]/2.0)*T
          - a[6]
        )*T
      + a[5]
    );
}


inline Foam::scalar Foam::janaf::dCpdT(const scalar T) const
{
    const coeffArray& a = coeffs(T);
    return ((4*a[4]*T + 3*a[3])*T + 2*a[2])*T + a[1];
}


inline Foam::scalar Foam::janaf::deltaS(const scalar T) const
{
    //TODO: Include contribution model
    const coeffArray& a = coeffs(T);
    return
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
      + a[6];
}


// ************************************************************************* //
