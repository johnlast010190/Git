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
    (c) 2016-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeChemistryReductionMethods_H
#define makeChemistryReductionMethods_H

#include "chemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethod.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistryReductionMethod(nullArg, ThermoPhysics)                 \
                                                                               \
    typedef chemistryReductionMethod<ThermoPhysics>                            \
        chemistryReductionMethod##ThermoPhysics;                               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethod##ThermoPhysics,                               \
        (                                                                      \
            word(chemistryReductionMethod##ThermoPhysics::typeName_())         \
          + "<" + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryReductionMethod##ThermoPhysics,                               \
        dictionary                                                             \
    )


#define makeChemistryReductionMethod(Method, ThermoPhysics)                    \
                                                                               \
    typedef chemistryReductionMethod<ThermoPhysics>                            \
        chemistryReductionMethod##ThermoPhysics;                               \
                                                                               \
    typedef chemistryReductionMethods::Method<ThermoPhysics>                   \
        chemistryReductionMethod##Method##ThermoPhysics;                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethod##Method##ThermoPhysics,                       \
        (                                                                      \
            word(chemistryReductionMethod##Method##ThermoPhysics::typeName_()) \
         + "<" + ThermoPhysics::typeName() + ">"                               \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        chemistryReductionMethod##ThermoPhysics,                               \
        chemistryReductionMethod##Method##ThermoPhysics,                       \
        dictionary                                                             \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
