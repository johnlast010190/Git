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

inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nCoeff_> coeffs
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    coeffs_(coeffs)
{}


inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta")),
    coeffs_(dict.lookup("coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::powerSeriesReactionRate::preEvaluate() const
{}


inline void Foam::powerSeriesReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::powerSeriesReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;

    forAll(coeffs_, n)
    {
        expArg += coeffs_[n]/pow(T, n + 1);
    }

    lta *= exp(expArg);

    return lta;
}


inline Foam::scalar Foam::powerSeriesReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0;
    scalar deriv = 0;

    forAll(coeffs_, n)
    {
        scalar cT = coeffs_[n]/pow(T, n + 1);
        expArg += cT;
        deriv -= (n + 1)*cT;
    }

    lta *= exp(expArg);

    return lta*(beta_ + deriv)/T;
}


inline bool Foam::powerSeriesReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::powerSeriesReactionRate::ddc
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


inline void Foam::powerSeriesReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("beta", beta_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("coeffs", coeffs_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerSeriesReactionRate& psrr
)
{
    psrr.write(os);
    return os;
}


// ************************************************************************* //
