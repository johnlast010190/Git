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
    (c) 2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::ArdenBuck::zeroC(273.15);
const Foam::scalar Foam::ArdenBuck::A(611.21);
const Foam::scalar Foam::ArdenBuck::B(18.678);
const Foam::scalar Foam::ArdenBuck::C(234.5);
const Foam::scalar Foam::ArdenBuck::D(257.14);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ArdenBuck::xByTC(scalar Tc) const
{
    return (B - Tc/C)/(D + Tc);
}


inline Foam::scalar Foam::ArdenBuck::pSat(const scalar T) const
{
    const scalar Tc = T - zeroC;

    return A*exp(Tc*xByTC(Tc));
}


inline Foam::scalar Foam::ArdenBuck::pSatPrime(const scalar T) const
{
    const scalar Tc = T - zeroC;
    const scalar x(xByTC(Tc));

    return A*exp(Tc*x)*(D*x - Tc/C)/(D + Tc);
}


inline Foam::scalar Foam::ArdenBuck::lnPSat(const scalar T) const
{
    const scalar Tc(T - zeroC);

    return log(A) + Tc*xByTC(Tc);
}


// ************************************************************************* //
