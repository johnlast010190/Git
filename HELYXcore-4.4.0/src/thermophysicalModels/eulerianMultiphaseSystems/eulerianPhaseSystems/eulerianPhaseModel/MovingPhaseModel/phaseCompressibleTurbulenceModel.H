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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Typedef
    Foam::phaseCompressibleTurbulenceModel

Description
    Typedef for phaseCompressibleTurbulenceModel

\*---------------------------------------------------------------------------*/

#ifndef phaseCompressibleTurbulenceModel_H
#define phaseCompressibleTurbulenceModel_H

#include "phaseCompressibleTurbulenceModelFwd.H"
#include "PhaseCompressibleTurbulenceModel/PhaseCompressibleTurbulenceModel.H"
#include "ThermalDiffusivity/ThermalDiffusivity.H"
#include "../eulerianPhaseModel/eulerianPhaseModel.H"

#include "turbulentFluidThermoModels/makeTurbulenceModel.H"

#include "laminar/laminarModel/laminarModel.H"
#include "RAS/RASModel/RASModel.H"
#include "LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
      typedef ThermalDiffusivity<PhaseCompressibleTurbulenceModel<eulerianPhaseModel>>
          phaseCompressibleTurbulenceModel;
}

makeTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    eulerianPhaseModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (eulerianPhaseModelPhaseCompressibleTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (eulerianPhaseModelPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (eulerianPhaseModelPhaseCompressibleTurbulenceModel, LES, Type)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
