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
    (c) 2023 Engys Ltd

\*---------------------------------------------------------------------------*/

inline const Foam::phasesTable&
Foam::phaseVolumeFractions::phases() const
{
    return phases_;
}


inline bool Foam::phaseVolumeFractions::contains
(
    const word& phaseName
) const
{
    return phases_.found(phaseName);
}


inline bool Foam::phaseVolumeFractions::active(label phasei) const
{
    return active_[phasei];
}


inline const Foam::List<bool>& Foam::phaseVolumeFractions::active() const
{
    return active_;
}


inline void Foam::phaseVolumeFractions::setActive(label phasei) const
{
    active_[phasei] = true;
    const_cast<volScalarField&>(alphas_[phasei]).writeOpt() =
        IOobject::AUTO_WRITE;
}


inline void Foam::phaseVolumeFractions::setInactive(label phasei) const
{
    active_[phasei] = false;
    const_cast<volScalarField&>(alphas_[phasei]).writeOpt() =
        IOobject::NO_WRITE;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::phaseVolumeFractions::alphas()
{
    return alphas_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::phaseVolumeFractions::alphas() const
{
    return alphas_;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::phaseVolumeFractions::fractions()
{
    return alphas();
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::phaseVolumeFractions::fractions() const
{
    return alphas();
}



inline Foam::volScalarField& Foam::phaseVolumeFractions::alphas(const label i)
{
    return alphas_[i];
}


inline const Foam::volScalarField& Foam::phaseVolumeFractions::alphas
(
    const label i
) const
{
    return alphas_[i];
}


inline Foam::volScalarField& Foam::phaseVolumeFractions::alphas
(
    const word& phaseName
)
{
    return alphas_[phases_[phaseName]];
}


inline const Foam::volScalarField& Foam::phaseVolumeFractions::alphas
(
    const word& phaseName
) const
{
    return alphas_[phases_[phaseName]];
}


// ************************************************************************* //
