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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef makeMulticomponentThermo_H
#define makeMulticomponentThermo_H

#include "mixtures/SpecieMixture/SpecieMixture.H"
#include "include/makeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineMulticomponentThermo(                                            \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)                    \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseMulticomponentThermo::composite,                                   \
        SpecieMixture<Mixture<ThermoPhys>>                                     \
    > CThermo##Mixture##ThermoPhys;                                            \
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


#define makeMulticomponentThermo(                                              \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)                    \
    defineMulticomponentThermo(                                                \
        BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys);               \
    addThermo(fluidMulticomponentThermo, CThermo, Mixture, ThermoPhys);        \
    addThermo(BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)


#define makeMulticomponentThermos(                                             \
    BaseThermo, BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)        \
    defineMulticomponentThermo(                                                \
        BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys);               \
    addThermo(basicThermo, CThermo, Mixture, ThermoPhys);                      \
    addThermo(fluidThermo, CThermo, Mixture, ThermoPhys);                      \
    addThermo(BaseThermo, CThermo, Mixture, ThermoPhys);                       \
    addThermo(fluidMulticomponentThermo, CThermo, Mixture, ThermoPhys);        \
    addThermo(BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys)


#define defineMaterialPhysicsMulticomponentThermo(                             \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys, MatType)           \
    typedef CThermo                                                            \
    <                                                                          \
        BaseMulticomponentThermo::composite,                                   \
        SpecieMixture<Mixture<ThermoPhys>>                                     \
    > MatType;                                                                 \
    defineTemplateTypeNameAndDebugWithName(MatType, #MatType, 0);


#define makeMaterialPhysicsMulticomponentThermos(                              \
    BaseThermo, BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys, MatType)\
    defineMaterialPhysicsMulticomponentThermo                                  \
    (                                                                          \
        BaseMulticomponentThermo,                                              \
        CThermo,                                                               \
        Mixture,                                                               \
        ThermoPhys,                                                            \
        MatType                                                                \
    )                                                                          \
    addToRunTimeSelectionTable(basicThermo, MatType, objectRegistry);          \
    addToRunTimeSelectionTable(fluidThermo, MatType, objectRegistry);          \
    addToRunTimeSelectionTable(BaseThermo, MatType, objectRegistry);           \
    addToRunTimeSelectionTable(                                                \
        fluidMulticomponentThermo, MatType, objectRegistry);                   \
    addToRunTimeSelectionTable(                                                \
        BaseMulticomponentThermo, MatType, objectRegistry);


#define makePremixedMaterials(Mixture, ThermoPhys, MatType)                    \
                                                                               \
    defineMaterialPhysicsMulticomponentThermo                                  \
    (                                                                          \
        psiuMulticomponentThermo,                                              \
        PsiuMulticomponentThermo,                                              \
        Mixture,                                                               \
        ThermoPhys,                                                            \
        MatType                                                                \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable(basicThermo, MatType, objectRegistry);          \
    addToRunTimeSelectionTable(fluidThermo, MatType, objectRegistry);          \
    addToRunTimeSelectionTable(rhoThermo, MatType, objectRegistry);            \
    addToRunTimeSelectionTable(psiuMulticomponentThermo, MatType, objectRegistry)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
