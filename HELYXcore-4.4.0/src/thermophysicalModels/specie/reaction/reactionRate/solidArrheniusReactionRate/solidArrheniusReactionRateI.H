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
    (c) 2011-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const scalar A,
    const scalar Ta,
    const scalar Tcrit
    //const scalar nReact
)
:
    A_(A),
    Ta_(Ta),
    Tcrit_(Tcrit)
{}


inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    Ta_(dict.lookup<scalar>("Ta")),
    Tcrit_(dict.lookup<scalar>("Tcrit"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::solidArrheniusReactionRate::preEvaluate() const
{}


inline void Foam::solidArrheniusReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::solidArrheniusReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField&,
    const label li
) const
{
    return (T < Tcrit_) ? 0 : A_*exp(-Ta_/T);
}


inline Foam::scalar Foam::solidArrheniusReactionRate::ddT
(
    const scalar,
    const scalar T,
    const scalarField&,
    const label
) const
{
    return (T < Tcrit_) ? 0 : A_*exp(-Ta_/T)*Ta_/sqr(T);
}


inline bool Foam::solidArrheniusReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::solidArrheniusReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    ddc = 0;
}


inline void Foam::solidArrheniusReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("Tcrit", Tcrit_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
