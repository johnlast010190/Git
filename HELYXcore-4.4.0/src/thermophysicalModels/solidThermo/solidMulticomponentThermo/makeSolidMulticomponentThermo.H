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
    Macros for creating reacting solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeReactingSolidThermo_H
#define makeReactingSolidThermo_H

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "mixtures/SpecieMixture/SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeSolidReactionMixtureThermo(BaseThermo,Cthermo,Mixture,ThermoPhys)  \
                                                                               \
    typedef SpecieMixture<Mixture<ThermoPhys>>                                 \
        SpecieMixture##Mixture##ThermoPhys;                                    \
                                                                               \
    typedef legacyBasicThermo                                                  \
    <                                                                          \
        BaseThermo::composite,                                                 \
        SpecieMixture##Mixture##ThermoPhys                                     \
    > legacyBasicThermo##Mixture##ThermoPhys;                                  \
                                                                               \
    typedef Cthermo                                                            \
    <                                                                          \
        BaseThermo::composite,                                                 \
        SpecieMixture##Mixture##ThermoPhys                                     \
    > Cthermo##Mixture##ThermoPhys;                                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Cthermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            #Cthermo"<" + Mixture<ThermoPhys>::typeName() + ">"                \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicThermo,                                                           \
        Cthermo##Mixture##ThermoPhys,                                          \
        objectRegistry                                                         \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        Cthermo##Mixture##ThermoPhys,                                          \
        objectRegistry                                                         \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        Cthermo##Mixture##ThermoPhys,                                          \
        dictionary                                                             \
    );


#define makeReactingSolidThermo(\
    BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)        \
                                                                               \
    typedef Transport                                                          \
    <                                                                          \
        species::thermo                                                        \
        <                                                                      \
            Thermo<EqnOfState<Specie>>,                                        \
            Type                                                               \
        >                                                                      \
    > Transport##Type##Thermo##EqnOfState##Specie;                             \
                                                                               \
    makeSolidReactionMixtureThermo                                             \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    );


#define makeMaterialPhysicsSolidMulticomponentThermos(                         \
    BaseMulticomponentThermo, CThermo, Mixture, ThermoPhys, MatType)           \
    typedef CThermo                                                            \
    <                                                                          \
        BaseMulticomponentThermo::composite,                                   \
        SpecieMixture<Mixture<ThermoPhys>>                                     \
    > solidMatType;                                                            \
    defineTemplateTypeNameAndDebugWithName(solidMatType, #MatType, 0);         \
    addToRunTimeSelectionTable(basicThermo, solidMatType, objectRegistry);     \
    addToRunTimeSelectionTable(BaseMulticomponentThermo, solidMatType, objectRegistry);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
