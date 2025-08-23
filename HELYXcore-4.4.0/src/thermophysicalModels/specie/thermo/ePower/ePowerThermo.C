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
    (c) 2020-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermo/ePower/ePowerThermo.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::ePowerThermo<EquationOfState>::ePowerThermo
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    EquationOfState(obr, dict),
    c0_(dict.subDict("thermodynamics").lookup<scalar>("C0")),
    n0_(dict.subDict("thermodynamics").lookup<scalar>("n0")),
    Tref_(dict.subDict("thermodynamics").lookup<scalar>("Tref")),
    hf_(dict.subDict("thermodynamics").lookup<scalar>("Hf"))
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ePowerThermo<EquationOfState>& et
)
{
    et.write(os);
    return os;
}


// ************************************************************************* //