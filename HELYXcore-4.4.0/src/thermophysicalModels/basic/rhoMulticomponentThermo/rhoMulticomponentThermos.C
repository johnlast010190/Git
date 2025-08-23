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
    (c) 2012-2023 OpenFOAM Foundation
    (c) 2017-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rhoThermo/rhoThermo.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "rhoThermo/legacyRhoThermo.H"
#include "mixtures/coefficientMulticomponentMixture/coefficientMulticomponentMixture.H"
#include "mixtures/valueMulticomponentMixture/valueMulticomponentMixture.H"
#include "rhoThermo/RhoThermo.H"
#include "include/forGases.H"
#include "include/forLiquids.H"
#include "include/forPolynomials.H"
#include "include/forTabulated.H"
#include "makeMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeRhoMulticomponentThermos(Mixture, ThermoPhysics)                   \
    makeMulticomponentThermos                                                  \
    (                                                                          \
        rhoThermo,                                                             \
        rhoMulticomponentThermo,                                               \
        legacyRhoThermo,                                                       \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

#define makeRhoMulticomponentThermo(Mixture, ThermoPhysics)                    \
    makeMulticomponentThermo                                                   \
    (                                                                          \
        rhoMulticomponentThermo,                                               \
        legacyRhoThermo,                                                       \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

namespace Foam
{
    forGases(makeRhoMulticomponentThermos, coefficientMulticomponentMixture);
    forLiquids(makeRhoMulticomponentThermos, coefficientMulticomponentMixture);
    forPolynomials(makeRhoMulticomponentThermos, coefficientMulticomponentMixture);
    forTabulated(makeRhoMulticomponentThermos, valueMulticomponentMixture);

    // New material library thermo
    makeMaterialPhysicsMulticomponentThermos
    (
        rhoThermo,
        rhoMulticomponentThermo,
        RhoThermo,
        multicomponentMixture,
        dummyThermo,
        reactingFluid
    );
}

// ************************************************************************* //
