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

inline Foam::thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const thirdBodyEfficiencies& tbes
)
:
    ArrheniusReactionRate(A, beta, Ta),
    thirdBodyEfficiencies_(tbes)
{}


inline Foam::thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dict),
    thirdBodyEfficiencies_(species, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::thirdBodyArrheniusReactionRate::preEvaluate() const
{
    ArrheniusReactionRate::preEvaluate();
}


inline void Foam::thirdBodyArrheniusReactionRate::postEvaluate() const
{
    ArrheniusReactionRate::postEvaluate();
}


inline Foam::scalar Foam::thirdBodyArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *ArrheniusReactionRate::operator()(p, T, c, li);
}


inline Foam::scalar Foam::thirdBodyArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *ArrheniusReactionRate::ddT(p, T, c, li);
}


inline bool Foam::thirdBodyArrheniusReactionRate::hasDdc() const
{
    return true;
}


inline void Foam::thirdBodyArrheniusReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    const scalar k = ArrheniusReactionRate::operator()(p, T, c, li);

    ddc = thirdBodyEfficiencies_.dMdc(c);
    ddc *= k;
}


inline void Foam::thirdBodyArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    thirdBodyEfficiencies_.write(os);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thirdBodyArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
