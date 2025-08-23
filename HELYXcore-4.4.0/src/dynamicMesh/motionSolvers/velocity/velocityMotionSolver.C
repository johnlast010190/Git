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

\*---------------------------------------------------------------------------*/

#include "motionSolvers/velocity/velocityMotionSolver.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityMotionSolver, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityMotionSolver::velocityMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityMotionSolver::~velocityMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityMotionSolver::movePoints(const pointField& p)
{
    // No local data that needs adapting.
}


void Foam::velocityMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    // fvMesh updates pointFields.
}


void Foam::velocityMotionSolver::mapMesh(const polyMeshMap& map)
{
    pointMotionU_.forceAssign(vector::zero);
    pointMotionU_.correctBoundaryConditions();
}


void Foam::velocityMotionSolver::distribute(const polyDistributionMap&)
{}


// ************************************************************************* //
