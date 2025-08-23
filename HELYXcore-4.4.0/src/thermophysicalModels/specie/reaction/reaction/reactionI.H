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
    (c) 2020 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "reaction.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::reaction::name() const
{
    return name_;
}


inline bool Foam::reaction::hasGaseousSpecies() const
{
    return gaseousSpcecies_;
}


inline const Foam::List<Foam::specieCoeffs>& Foam::reaction::lhs() const
{
    return lhs_;
}


inline const Foam::List<Foam::specieCoeffs>& Foam::reaction::rhs() const
{
    return rhs_;
}


const Foam::List<Foam::specieCoeffs>& Foam::reaction::glhs() const
{
    return glhs_;
}


const Foam::List<Foam::specieCoeffs>& Foam::reaction::grhs() const
{
    return grhs_;
}


const Foam::speciesTable& Foam::reaction::gasSpecies() const
{
    return pyrolisisGases_;
}


inline const Foam::speciesTable& Foam::reaction::species() const
{
    return species_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const reaction& r)
{
    r.write(os);
    return os;
}


// ************************************************************************* //