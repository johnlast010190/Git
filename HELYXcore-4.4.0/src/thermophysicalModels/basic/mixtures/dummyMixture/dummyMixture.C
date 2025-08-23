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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "mixtures/dummyMixture/dummyMixture.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::dummyMixture<ThermoType>::dummyMixture
(
    const dictionary& thermoDict,
    const objectRegistry& mesh,
    const word& phaseName
)
:
    basicCombustionMixture
    (
        thermoDict,
        speciesTable(thermoDict.lookup<wordList>("species")),
        mesh,
        phaseName
    ),
    stoicRatio_("stoichiometricAirFuelMassRatio", dimless, thermoDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::dummyMixture<ThermoType>::read
(
    const objectRegistry& mesh,
    const dictionary& thermoDict
)
{
    thermoDict.lookup("stoichiometricAirFuelMassRatio") >> stoicRatio_;
}


template<class ThermoType>
const ThermoType& Foam::dummyMixture<ThermoType>::specieThermo
(
    const label speciei
) const
{
    NotImplemented;
}


// ************************************************************************* //
