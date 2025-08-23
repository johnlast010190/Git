
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
    (c) 2020-2021 OpenFOAM Foundation

InClass
    Foam::fluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "fluidThermo/fluidThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermo(BaseThermo, CThermo, Mixture, ThermoPhys)                 \
                                                                               \
    typedef CThermo<BaseThermo::composite, Mixture<ThermoPhys>>                \
        CThermo##Mixture##ThermoPhys;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            word(CThermo##Mixture##ThermoPhys::typeName_()) + "<"              \
          + Mixture<ThermoPhys>::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define addThermo(BaseThermo, CThermo, Mixture, ThermoPhys)                    \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo##Mixture##ThermoPhys,                                          \
        objectRegistry                                                         \
    )


#define makeThermo(BaseThermo, CThermo, Mixture, ThermoPhys)                   \
    defineThermo(BaseThermo, CThermo, Mixture, ThermoPhys);                    \
    addThermo(basicThermo, CThermo, Mixture, ThermoPhys);                      \
    addThermo(fluidThermo, CThermo, Mixture, ThermoPhys);                      \
    addThermo(BaseThermo, CThermo, Mixture, ThermoPhys)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //