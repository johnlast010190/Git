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

#ifndef forGases_H
#define forGases_H

#include "specie/specie.H"
#include "equationOfState/Boussinesq/Boussinesq.H"
#include "equationOfState/incompressiblePerfectGas/incompressiblePerfectGas.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "equationOfState/PengRobinsonGas/PengRobinsonGas.H"
#include "thermo/eConst/eConstThermo.H"
#include "thermo/hConst/hConstThermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "transport/const/constTransport.H"
#include "transport/sutherland/sutherlandTransport.H"
#include "thermo/thermo/thermo.H"
#include "include/dummyThermo.H"
#include "include/forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forGasEquations(Mu, He, Cp, Macro, Args...)                            \
    forThermo(Mu, He, Cp, Boussinesq, specie, Macro, Args);                    \
    forThermo(Mu, He, Cp, incompressiblePerfectGas, specie, Macro, Args);      \
    forThermo(Mu, He, Cp, PengRobinsonGas, specie, Macro, Args);               \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forGasEnergiesAndThermos(Mu, Macro, Args...)                           \
    forGasEquations(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);          \
    forGasEquations(Mu, sensibleEnthalpy, janafThermo, Macro, Args);           \
    forGasEquations(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);    \
    forGasEquations(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args);    \
    forGasEquations(Mu, sensibleInternalEnergy, janafThermo, Macro, Args)

#define forGasTransports(Macro, Args...)                                       \
    forGasEnergiesAndThermos(constTransport, Macro, Args);                     \
    forGasEnergiesAndThermos(sutherlandTransport, Macro, Args)

#define forGases(Macro, Args...)                                               \
    forGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //