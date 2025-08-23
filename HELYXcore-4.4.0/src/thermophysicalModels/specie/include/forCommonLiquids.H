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

#ifndef forCommonLiquids_H
#define forCommonLiquids_H

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

#define forCommonLiquidEqns(Mu, He, Cp, Macro, Args...)                        \
    forThermo(Mu, He, Cp, rhoConst, specie, Macro, Args);

#define forCommonLiquidEnergiesAndThermos(Mu, Macro, Args...)                  \
    forCommonLiquidEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);      \
    forCommonLiquidEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);\
    forCommonLiquidEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forCommonLiquidTransports(Macro, Args...)                              \
    forCommonLiquidEnergiesAndThermos(constTransport, Macro, Args);

#define forCommonLiquids(Macro, Args...)                                       \
    forCommonLiquidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //