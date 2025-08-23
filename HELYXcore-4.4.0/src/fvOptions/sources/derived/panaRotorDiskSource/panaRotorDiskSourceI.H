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
    (c) 2012-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "panaRotorDiskSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::fv::panaRotorDiskSource::rhoRef() const
{
    return rhoRef_;
}


Foam::scalar Foam::fv::panaRotorDiskSource::omega() const
{
    return omega_;
}


const Foam::List<Foam::point>& Foam::fv::panaRotorDiskSource::x() const
{
    return x_;
}


const Foam::cylindricalCS& Foam::fv::panaRotorDiskSource::coordSys() const
{
    return coordSys_;
}


bool Foam::fv::panaRotorDiskSource::compressible() const
{
    return rhoName_ != "none";
}


Foam::tmp<Foam::volScalarField> Foam::fv::panaRotorDiskSource::rho() const
{
    if (compressible())
    {
        return obr_.lookupObject<volScalarField>(rhoName_);
    }
    else
    {
        return volScalarField::null();
    }
}


// ************************************************************************* //
