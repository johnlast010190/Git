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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef makeChemistrySolver_H
#define makeChemistrySolver_H

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolver(Model, ThermoPhysics)                            \
    typedef Model<ThermoPhysics> Model##ThermoPhysics;                         \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ThermoPhysics,                                                  \
        (                                                                      \
            word(Model##ThermoPhysics::typeName_())                            \
          + "<" + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    );


#define makeChemistrySolver(Solver, Model, ThermoPhysics)                      \
                                                                               \
    typedef Solver<Model<ThermoPhysics>> Solver##Model##ThermoPhysics;         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##ThermoPhysics,                                          \
        (                                                                      \
            word(Solver##Model##ThermoPhysics::typeName_())                    \
          + "<" + word(Model<ThermoPhysics>::typeName_())                      \
          + "<" + ThermoPhysics::typeName() + ">>"                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicChemistryModel,                                                   \
        Solver##Model##ThermoPhysics,                                          \
        thermo                                                                 \
    )


#define defineChemistrySolvers(nullArg, ThermoPhysics)                         \
    defineChemistrySolver(chemistryModel, ThermoPhysics)


#define makeChemistrySolvers(Solver, ThermoPhysics)                            \
    makeChemistrySolver(Solver, chemistryModel, ThermoPhysics)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
