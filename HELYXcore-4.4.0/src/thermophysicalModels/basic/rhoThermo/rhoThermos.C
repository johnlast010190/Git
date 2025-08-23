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
    (c) 2017-2024 Engys Ltd.
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "rhoThermo/rhoThermo.H"
#include "rhoThermo/legacyRhoThermo.H"
#include "mixtures/pureMixture/pureMixture.H"

#include "include/forGases.H"
#include "include/forLiquids.H"
#include "include/forPolynomials.H"
#include "include/makeThermo.H"
#include "include/forTabulated.H"

#include "materialModels/materialMacros.H"
#include "rhoThermo/RhoThermo.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forGases(makeThermo, rhoThermo, legacyRhoThermo, pureMixture);
    forLiquids(makeThermo, rhoThermo, legacyRhoThermo, pureMixture);
    forPolynomials(makeThermo, rhoThermo, legacyRhoThermo, pureMixture);
    forTabulated(makeThermo, rhoThermo, legacyRhoThermo, pureMixture);
    makeMatThermo(rhoThermo, speciesMassFractions, RhoThermo, fluid);
}

// ************************************************************************* //
