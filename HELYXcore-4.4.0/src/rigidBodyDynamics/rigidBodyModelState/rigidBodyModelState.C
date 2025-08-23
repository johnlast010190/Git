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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "rigidBodyModelState/rigidBodyModelState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyModelState::rigidBodyModelState
(
    const rigidBodyModel& model
)
:
    q_(model.nDoF(), Zero),
    qDot_(model.nDoF(), Zero),
    qDdot_(model.nDoF(), Zero),
    deltaT_(0)
{}


Foam::RBD::rigidBodyModelState::rigidBodyModelState
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    q_(dict.lookupOrDefault("q", scalarField(model.nDoF(), Zero))),
    qDot_(dict.lookupOrDefault("qDot", scalarField(model.nDoF(), Zero))),
    qDdot_(dict.lookupOrDefault("qDdot", scalarField(model.nDoF(), Zero))),
    deltaT_(dict.lookupOrDefault<scalar>("deltaT", 0))
{}


// ************************************************************************* //
