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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicMotionSolverFvMesh/dynamicMotionSolverFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/volFields/volFields.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverFvMesh::dynamicMotionSolverFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    motionPtr_(motionSolver::New(*this, dynamicMeshDict()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverFvMesh::~dynamicMotionSolverFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::dynamicMotionSolverFvMesh::motion() const
{
    return motionPtr_();
}


bool Foam::dynamicMotionSolverFvMesh::update()
{
    // Update frame motions
    coordinateFrame::updateStates((*this).thisDb());

    fvMesh::movePoints(motionPtr_->newPoints());

    if (foundObject<volVectorField>("U"))
    {
        volVectorField& U =
            const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
        U.correctBoundaryConditions();
    }

    return true;
}


bool Foam::dynamicMotionSolverFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    motionPtr_->write();

    return fvMesh::writeObject(fmt, ver, cmp, write);
}


// ************************************************************************* //
