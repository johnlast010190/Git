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

inline Foam::JanevReactionRate::JanevReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nb_> b
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    b_(b)
{}


inline Foam::JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta")),
    b_(dict.lookup("b"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::JanevReactionRate::preEvaluate() const
{}


inline void Foam::JanevReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::JanevReactionRate::operator()
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

    if (mag(Ta_) > VSMALL)
    {
        expArg -= Ta_/T;
    }

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    lta *= exp(expArg);

    return lta;
}


inline Foam::scalar Foam::JanevReactionRate::ddT
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

    if (mag(Ta_) > VSMALL)
    {
        expArg -= Ta_/T;
    }

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    scalar deriv = b_[1];

    for (int n=2; n<nb_; n++)
    {
        deriv += n*b_[n]*pow(lnT, n-1);
    }

    lta *= exp(expArg);

    return lta*(beta_ + Ta_/T + deriv)/T;
}


inline bool Foam::JanevReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::JanevReactionRate::ddc
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


inline void Foam::JanevReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << nl;
    os.writeKeyword("beta") << beta_ << nl;
    os.writeKeyword("Ta") << Ta_ << nl;
    os.writeKeyword("b") << b_ << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const JanevReactionRate& jrr
)
{
    jrr.write(os);
    return os;
}


// ************************************************************************* //
