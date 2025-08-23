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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "specie/specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidProperties::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::solidProperties::Cp() const
{
    return Cp_;
}


inline Foam::scalar Foam::solidProperties::kappa() const
{
    return kappa_;
}


inline Foam::scalar Foam::solidProperties::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::solidProperties::hs(const scalar T) const
{
    return Cp_*(T - Tstd);
}


inline Foam::scalar Foam::solidProperties::ha(const scalar T) const
{
    return hs(T) + hf();
}


inline Foam::scalar Foam::solidProperties::emissivity() const
{
    return emissivity_;
}


// ************************************************************************* //
