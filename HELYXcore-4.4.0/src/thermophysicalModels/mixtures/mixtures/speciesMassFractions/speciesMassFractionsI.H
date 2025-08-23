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
    (c) 2011-2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

inline const Foam::speciesTable&
Foam::speciesMassFractions::species() const
{
    return species_;
}


inline bool Foam::speciesMassFractions::contains(const word& specieName) const
{
    return species_.found(specieName);
}


inline Foam::label Foam::speciesMassFractions::defaultSpecie() const
{
    return defaultSpecieIndex_;
}


inline bool Foam::speciesMassFractions::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::speciesMassFractions::active() const
{
    return active_;
}


inline void Foam::speciesMassFractions::setActive(label speciei) const
{
    active_[speciei] = true;
    const_cast<volScalarField&>(Y_[speciei]).writeOpt() = IOobject::AUTO_WRITE;
}


inline void Foam::speciesMassFractions::setInactive(label speciei) const
{
    active_[speciei] = false;
    const_cast<volScalarField&>(Y_[speciei]).writeOpt() = IOobject::NO_WRITE;
}


inline bool Foam::speciesMassFractions::solve(label speciei) const
{
    return speciei != defaultSpecieIndex_ && active_[speciei];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::speciesMassFractions::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::speciesMassFractions::Y() const
{
    return Y_;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::speciesMassFractions::fractions()
{
    return Y();
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::speciesMassFractions::fractions() const
{
    return Y();
}



inline Foam::volScalarField& Foam::speciesMassFractions::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::speciesMassFractions::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::speciesMassFractions::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::speciesMassFractions::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


inline Foam::label Foam::speciesMassFractions::index
(
    const volScalarField& Yi
) const
{
    return species_[Yi.member()];
}


// ************************************************************************* //
