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
    (c) 2018-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::MichaelisMentenReactionRate::MichaelisMentenReactionRate
(
    const speciesTable& st,
    const dictionary& dict
)
:
    species_(st),
    Vmax_(readScalar(dict.lookup("Vmax"))),
    Km_(readScalar(dict.lookup("Km"))),
    s_(st[dict.lookup("S")])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::MichaelisMentenReactionRate::preEvaluate() const
{}


inline void Foam::MichaelisMentenReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::MichaelisMentenReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return Vmax_/(Km_ + c[s_]);
}


inline Foam::scalar Foam::MichaelisMentenReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline bool Foam::MichaelisMentenReactionRate::hasDdc() const
{
    return true;
}


inline void Foam::MichaelisMentenReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    ddc = 0;
    ddc[s_] = - Vmax_/sqr(Km_ + c[s_]);
}


inline void Foam::MichaelisMentenReactionRate::write(Ostream& os) const
{
    os.writeKeyword("Vmax") << Vmax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Km") << Km_ << token::END_STATEMENT << nl;
    os.writeKeyword("S") << species_[s_] << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MichaelisMentenReactionRate& mmrr
)
{
    mmrr.write(os);
    return os;
}


// ************************************************************************* //