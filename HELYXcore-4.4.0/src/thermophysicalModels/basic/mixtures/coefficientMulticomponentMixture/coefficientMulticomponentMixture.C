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
    (c) 2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/coefficientMulticomponentMixture/coefficientMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::coefficientMulticomponentMixture<ThermoType>::
coefficientMulticomponentMixture
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    multicomponentMixture<ThermoType>(thermoDict, obr, phaseName),
    mixture_("mixture", this->specieThermos()[0])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::cellThermoMixture
(
    const label celli
) const
{
    mixture_ = this->Y()[0][celli]*this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ += this->Y()[i][celli]*this->specieThermos()[i];
    }

    return mixture_;
}


template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::patchFaceThermoMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ =
        this->Y()[0].boundaryField()[patchi][facei]
       *this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ +=
            this->Y()[i].boundaryField()[patchi][facei]
           *this->specieThermos()[i];
    }

    return mixture_;
}


// ************************************************************************* //