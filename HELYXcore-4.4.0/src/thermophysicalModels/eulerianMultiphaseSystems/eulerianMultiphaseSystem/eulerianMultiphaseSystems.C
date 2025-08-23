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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "../eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "eulerianMultiphaseSystem.H"
#include "../eulerianPhaseSystems/PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.H"
#include "../eulerianPhaseSystems/PhaseSystems/HeatTransferPhaseSystem/HeatTransferPhaseSystem.H"
#include "../eulerianPhaseSystems/PhaseSystems/InterfaceCompositionPhaseChangePhaseSystem/InterfaceCompositionPhaseChangePhaseSystem.H"
#include "../eulerianPhaseSystems/PhaseSystems/ThermalPhaseChangePhaseSystem/ThermalPhaseChangePhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        HeatTransferPhaseSystem
        <
            MomentumTransferPhaseSystem<eulerianMultiphaseSystem>
        >
        heatAndMomentumTransferMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        eulerianMultiphaseSystem,
        heatAndMomentumTransferMultiphaseSystem,
        dictionary,
        heatAndMomentumTransfer
    );

    typedef
        InterfaceCompositionPhaseChangePhaseSystem
        <
            MomentumTransferPhaseSystem<eulerianMultiphaseSystem>
        >
        interfaceCompositionPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        eulerianMultiphaseSystem,
        interfaceCompositionPhaseChangeMultiphaseSystem,
        dictionary,
        interfaceCompositionPhaseChange
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            MomentumTransferPhaseSystem<eulerianMultiphaseSystem>
        >
        thermalPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        eulerianMultiphaseSystem,
        thermalPhaseChangeMultiphaseSystem,
        dictionary,
        thermalPhaseChange
    );
}


// ************************************************************************* //
