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

inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar B,
    const scalar C
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    B_(B),
    C_(C)
{}


inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta")),
    B_(dict.lookup<scalar>("B")),
    C_(dict.lookup<scalar>("C"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::LandauTellerReactionRate::preEvaluate() const
{}


inline void Foam::LandauTellerReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::LandauTellerReactionRate::operator()
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

    if (mag(B_) > VSMALL)
    {
        expArg += B_/cbrt(T);
    }

    if (mag(C_) > VSMALL)
    {
        expArg += C_/pow(T, 2.0/3.0);
    }

    if (mag(expArg) > VSMALL)
    {
        lta *= exp(expArg);
    }

    return lta;
}


inline Foam::scalar Foam::LandauTellerReactionRate::ddT
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

    if (mag(Ta_) > VSMALL)
    {
        scalar TaT = Ta_/T;
        expArg -= TaT;
        deriv += TaT;
    }

    if (mag(B_) > VSMALL)
    {
        scalar BT = B_/cbrt(T);
        expArg += BT;
        deriv -= BT/3;
    }

    if (mag(C_) > VSMALL)
    {
        scalar CT = C_/pow(T, 2.0/3.0);
        expArg += CT;
        deriv -= 2*CT/3;
    }

    if (mag(expArg) > VSMALL)
    {
        lta *= exp(expArg);
    }

    return lta*(beta_+deriv)/T;
}


inline bool Foam::LandauTellerReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::LandauTellerReactionRate::ddc
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


inline void Foam::LandauTellerReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("beta", beta_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("B", B_);
    os.writeEntry("C", C_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LandauTellerReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
