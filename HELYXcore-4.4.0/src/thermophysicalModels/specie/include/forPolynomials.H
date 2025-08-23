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

#ifndef forPolynomials_H
#define forPolynomials_H

#include "specie/specie.H"
#include "equationOfState/icoPolynomial/icoPolynomial.H"
#include "thermo/hPolynomial/hPolynomialThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "transport/logPolynomial/logPolynomialTransport.H"
#include "transport/polynomial/polynomialTransport.H"
#include "thermo/thermo/thermo.H"
#include "include/forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forPolynomialEquations(Mu, He, Cp, Args...)                            \
    forThermo(Mu, He, Cp, icoPolynomial, specie, Args);

#define forPolynomialEnergiesAndThermos(Mu, Args...)                           \
    forPolynomialEquations(Mu, sensibleEnthalpy, hPolynomialThermo, Args);     \
    forPolynomialEquations(Mu, sensibleInternalEnergy, hPolynomialThermo, Args)

#define forPolynomialTransports(Args...)                                       \
    forPolynomialEnergiesAndThermos(logPolynomialTransport, Args);             \
    forPolynomialEnergiesAndThermos(polynomialTransport, Args)

#define forPolynomials(Macro, Args...)                                         \
    forPolynomialTransports(Macro, Args)

#endif

// ************************************************************************* //