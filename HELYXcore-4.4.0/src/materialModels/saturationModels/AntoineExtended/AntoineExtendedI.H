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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::AntoineExtended::pSat(const scalar T) const
{
    return exp(A_ + B_/(C_ + T) + E_*pow(T, F_))*pow(T, D_);
}


inline Foam::scalar Foam::AntoineExtended::pSatPrime(const scalar T) const
{
    return pSat(T)*((D_ + E_*F_*pow(T, F_))/T - B_/sqr(C_ + T));
}


inline Foam::scalar Foam::AntoineExtended::lnPSat(const scalar T) const
{
    return A_ + B_/(C_ + T) + D_*log(T) + E_*pow(T, F_);
}


// ************************************************************************* //
