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

#ifndef forTabulated_H
#define forTabulated_H

#include "specie/specie.H"
#include "equationOfState/icoTabulated/icoTabulated.H"
#include "thermo/hTabulated/hTabulatedThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "transport/tabulated/tabulatedTransport.H"
#include "thermo/thermo/thermo.H"
#include "include/forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forTabulatedEqns(Mu, He, Cp, Macro, Args...)                          \
    forThermo(Mu, He, Cp, icoTabulated, specie, Macro, Args);

#define forTabulatedEnergiesAndThermos(Mu, Macro, Args...)                    \
    forTabulatedEqns(Mu, sensibleEnthalpy, hTabulatedThermo, Macro, Args);     \
    forTabulatedEqns(Mu, sensibleInternalEnergy, hTabulatedThermo, Macro, Args);

#define forTabulatedTransport(Macro, Args...)                                 \
    forTabulatedEnergiesAndThermos(tabulatedTransport, Macro, Args);

#define forTabulated(Macro, Args...)                                          \
    forTabulatedTransport(Macro, Args)

#endif

// ************************************************************************* //