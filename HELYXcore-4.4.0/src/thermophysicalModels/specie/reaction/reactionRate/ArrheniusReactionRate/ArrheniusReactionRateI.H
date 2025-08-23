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

inline Foam::ArrheniusReactionRate::ArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta)
{}


inline Foam::ArrheniusReactionRate::ArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::ArrheniusReactionRate::preEvaluate() const
{}


inline void Foam::ArrheniusReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::ArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label li
) const
{
    scalar ak = A_;

    if (mag(beta_) > VSMALL)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > VSMALL)
    {
        ak *= exp(-Ta_/T);
    }

    return ak;
}

inline Foam::scalar Foam::ArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label li
) const
{
    scalar ak = A_;

    if (mag(beta_) > VSMALL)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > VSMALL)
    {
        ak *= exp(-Ta_/T);
    }

    return ak*(beta_ + Ta_/T)/T;
}


inline bool Foam::ArrheniusReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::ArrheniusReactionRate::ddc
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


inline void Foam::ArrheniusReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("beta", beta_);
    os.writeEntry("Ta", Ta_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
