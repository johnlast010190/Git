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

#ifndef forLiquids_H
#define forLiquids_H

#include "specie/specie.H"
#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"
#include "equationOfState/rhoConst/rhoConst.H"
#include "thermo/eConst/eConstThermo.H"
#include "thermo/hConst/hConstThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "transport/const/constTransport.H"
#include "thermo/thermo/thermo.H"
#include "include/forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forLiquidEquations(Mu, He, Cp, Macro, Args...)                         \
    forThermo(Mu, He, Cp, adiabaticPerfectFluid, specie, Macro, Args);         \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args)

#define forLiquidEnergiesAndThermos(Mu, Macro, Args...)                        \
    forLiquidEquations(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);       \
    forLiquidEquations(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args); \
    forLiquidEquations(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forLiquidTransports(Macro, Args...)                                    \
    forLiquidEnergiesAndThermos(constTransport, Macro, Args)

#define forLiquids(Macro, Args...)                                             \
    forLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //