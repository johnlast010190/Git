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
    (c) 2015-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "rhoThermo/rhoThermo.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"

#include "combustionModel/combustionModel.H"

#include "eulerianPhaseModel.H"
#include "../ThermoPhaseModel/ThermoPhaseModel.H"
#include "../AnisothermalPhaseModel/AnisothermalPhaseModel.H"
#include "../PurePhaseModel/PurePhaseModel.H"
#include "../MulticomponentPhaseModel/MulticomponentPhaseModel.H"
#include "../InertPhaseModel/InertPhaseModel.H"
#include "../ReactingPhaseModel/ReactingPhaseModel.H"
#include "../MovingPhaseModel/MovingPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<eulerianPhaseModel, rhoThermo>
                    >
                >
            >
        >
        purePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        eulerianPhaseModel,
        purePhaseModel,
        eulerianPhaseSystem,
        pure
    );

    typedef
        AnisothermalPhaseModel
        <
            MulticomponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<eulerianPhaseModel, rhoMulticomponentThermo>
                    >
                >
            >
        >
        multicomponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        eulerianPhaseModel,
        multicomponentPhaseModel,
        eulerianPhaseSystem,
        multicomponent
    );

    typedef
        AnisothermalPhaseModel
        <
            MulticomponentPhaseModel
            <
                ReactingPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<eulerianPhaseModel, rhoMulticomponentThermo>
                    >
                >
            >
        >
        reactingPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        eulerianPhaseModel,
        reactingPhaseModel,
        eulerianPhaseSystem,
        reacting
    );
}

// ************************************************************************* //
