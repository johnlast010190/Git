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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hTabulated> Foam::hTabulated::clone() const
{
    return autoPtr<hTabulated>(new hTabulated(*this));
}


inline Foam::autoPtr<Foam::hTabulated> Foam::hTabulated::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hTabulated>
    (
        new hTabulated(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hTabulated::Cp
(
    const scalar p,
    const scalar T,
    const scalar CpContribution
) const
{
    return CpContribution + CpTable_()(p, T);
}


inline Foam::scalar Foam::hTabulated::ha
(
    const scalar p,
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + CpTable_().integrateOverY(p, 0, T) + hf_;
}


inline Foam::scalar Foam::hTabulated::hs
(
    const scalar p,
    const scalar T,
    const scalar hContribution
) const
{
    return hContribution + CpTable_().integrateOverY(p, 0, T);
}


inline Foam::scalar Foam::hTabulated::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::hTabulated::s
(
    const scalar p,
    const scalar T,
    const scalar spContribution
) const
{
    NotImplemented;
}


inline Foam::scalar Foam::hTabulated::dCpdT(const scalar T) const
{
    NotImplemented;
}


inline Foam::scalar Foam::hTabulated::Cv
(
    const scalar Cp,
    const scalar CpMCv
) const
{
    return Cp - CpMCv;
}


inline Foam::scalar Foam::hTabulated::es
(
    const scalar p,
    const scalar rho,
    const scalar hs
) const
{
    return hs - p/rho;
}


inline Foam::scalar Foam::hTabulated::ea
(
    const scalar p,
    const scalar rho,
    const scalar ha
) const
{
    return ha - p/rho;
}


// ************************************************************************* //
