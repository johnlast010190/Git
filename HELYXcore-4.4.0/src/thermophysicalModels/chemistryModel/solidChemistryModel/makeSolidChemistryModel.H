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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Description
    Macros for instantiating solid chemistry models

\*---------------------------------------------------------------------------*/

#ifndef makeSolidChemistryModel_H
#define makeSolidChemistryModel_H

#include "pyrolysisChemistryModel/pyrolysisChemistryModel.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "solidChemistryModel/solidChemistryModel.H"
#include "chemistrySolver/noChemistrySolver/noChemistrySolver.H"
#include "chemistrySolver/ode/ode.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

#define makeSolidGasChemistryModel(GThermo, SThermo)                           \
                                                                               \
    typedef solidChemistryModel<SThermo>                                       \
        solidChemistryModel##SThermo;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        solidChemistryModel##SThermo,                                          \
        (                                                                      \
            word(solidChemistryModel##SThermo::typeName_())                    \
          + "<" + SThermo::typeName() + ">"                                    \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
    typedef pyrolysisChemistryModel<SThermo, GThermo>                          \
        pyrolysisChemistryModel##SThermo##GThermo;                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        pyrolysisChemistryModel##SThermo##GThermo,                             \
        (                                                                      \
            word(pyrolysisChemistryModel##SThermo##GThermo::typeName_())       \
          + "<" + SThermo::typeName() + "," + GThermo::typeName() + ">"        \
        ).c_str(),                                                             \
        0                                                                      \
    );
}


#define makeSolidChemistrySolver(Solver, GThermo, SThermo)                     \
                                                                               \
    typedef Solver<pyrolysisChemistryModel<SThermo, GThermo>>                  \
        Solver##pyrolysisChemistryModel##SThermo##GThermo;                     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##pyrolysisChemistryModel##SThermo##GThermo,                     \
        (                                                                      \
            #Solver"<"                                                         \
          + word(pyrolysisChemistryModel<SThermo, GThermo>::typeName_())       \
          + "<" + SThermo::typeName() + "," + GThermo::typeName() + ">>"       \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicSolidChemistryModel,                                              \
        Solver##pyrolysisChemistryModel##SThermo##GThermo,                     \
        thermo                                                                 \
    );


#define makeSolidChemistrySolvers(GThermo, SThermo)                            \
    makeSolidChemistrySolver(noChemistrySolver, GThermo, SThermo);             \
    makeSolidChemistrySolver(ode, GThermo, SThermo);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
