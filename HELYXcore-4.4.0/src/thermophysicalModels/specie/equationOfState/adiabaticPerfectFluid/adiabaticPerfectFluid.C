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
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    Specie(obr, dict),
    p0_(dict.subDict("equationOfState").lookup<scalar>("p0")),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    gamma_(dict.subDict("equationOfState").lookup<scalar>("gamma")),
    B_(dict.subDict("equationOfState").lookup<scalar>("B"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::adiabaticPerfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("p0", p0_);
    dict.add("rho0", rho0_);
    dict.add("gamma", gamma_);
    dict.add("B", B_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const adiabaticPerfectFluid<Specie>& pf
)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
