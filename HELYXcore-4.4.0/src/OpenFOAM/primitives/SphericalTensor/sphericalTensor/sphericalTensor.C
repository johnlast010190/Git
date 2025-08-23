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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::sphericalTensor::vsType::typeName = "sphericalTensor";

template<>
const char* const Foam::sphericalTensor::vsType::componentNames[] = {"ii"};

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::zero
(
    sphericalTensor::uniform(0)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::one
(
    sphericalTensor::uniform(1)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::max
(
    sphericalTensor::uniform(VGREAT)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::min
(
    sphericalTensor::uniform(-VGREAT)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::rootMax
(
    sphericalTensor::uniform(ROOTVGREAT)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::rootMin
(
    sphericalTensor::uniform(-ROOTVGREAT)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::I(1);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::oneThirdI(1.0/3.0);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::twoThirdsI(2.0/3.0);

// ************************************************************************* //
