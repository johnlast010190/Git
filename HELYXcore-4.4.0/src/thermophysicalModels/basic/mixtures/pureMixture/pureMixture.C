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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/pureMixture/pureMixture.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::pureMixture<ThermoType>::pureMixture
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    basicMixture(thermoDict, obr, phaseName),
    mixture_(obr, thermoDict.subDict("mixture"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::pureMixture<ThermoType>::read
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    mixture_ = ThermoType(obr, thermoDict.subDict("mixture"));
}


// ************************************************************************* //
