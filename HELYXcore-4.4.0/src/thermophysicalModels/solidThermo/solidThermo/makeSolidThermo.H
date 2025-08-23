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

InClass
    Foam::solidThermo

Description
    Macros for creating solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeHeSolidThermo_H
#define makeHeSolidThermo_H

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys)            \
    typedef                                                                    \
        legacyBasicThermo                                                      \
        <                                                                      \
            BaseThermo::composite,                                             \
            Mixture<ThermoPhys>                                                \
        > legacyBasicThermo##Mixture##ThermoPhys;                              \
                                                                               \
    typedef                                                                    \
        CThermo                                                                \
        <                                                                      \
            BaseThermo::composite,                                             \
            Mixture<ThermoPhys>                                                \
        > CThermo##Mixture##ThermoPhys;                                        \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            #CThermo"<"#Mixture"<"                                             \
          + ThermoPhys::typeName()                                             \
          + ">>"                                                               \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define addSolidThermo(BaseThermo, CThermoMixtureThermoPhys, Table)            \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        Table                                                                  \
    )


#define makeSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys)              \
    defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
    addSolidThermo(BaseThermo, CThermo##Mixture##ThermoPhys, objectRegistry)


#define makeSolidThermos(BaseThermo, CThermo, Mixture, ThermoPhys)             \
    defineSolidThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
    addSolidThermo(basicThermo, CThermo##Mixture##ThermoPhys, objectRegistry); \
    addSolidThermo(BaseThermo, CThermo##Mixture##ThermoPhys, objectRegistry);  \
    addSolidThermo(BaseThermo, CThermo##Mixture##ThermoPhys, dictionary)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
