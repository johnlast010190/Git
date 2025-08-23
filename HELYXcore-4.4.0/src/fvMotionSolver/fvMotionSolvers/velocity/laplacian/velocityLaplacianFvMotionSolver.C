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
    (c) 2013-2024 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMotionSolvers/velocity/laplacian/velocityLaplacianFvMotionSolver.H"
#include "motionInterpolation/motionInterpolation/motionInterpolation.H"
#include "motionDiffusivity/motionDiffusivity/motionDiffusivity.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMatrices/fvMatrices.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "motionInterpolation/patchCorrected/patchCorrectedInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        velocityLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityLaplacianFvMotionSolver::velocityLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    velocityMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
     ),
    fixedVelocity_(mesh.nPoints(), vector(GREAT, GREAT, GREAT))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityLaplacianFvMotionSolver::~velocityLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityLaplacianFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate(cellMotionU_, pointMotionU_);

    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + fvMesh_.time().deltaTValue()*pointMotionU_.primitiveField()
    );

    pointField& cPts = tcurPoints.ref();

    // The processor boundaries go out of sync due to the interpolation
    // hence need to sync them. Should be handled by more consitent
    // interpolation across the processor boundaries
    if
    (
        Pstream::parRun()
     && interpolationPtr_->type() == patchCorrectedInterpolation::typeName
    )
    {
        pointField p1(cPts);
        pointField p2(cPts);
        const vector greatPoint(GREAT, GREAT, GREAT);
        syncTools::syncPointPositions
        (
            fvMesh_,
            p1,
            maxEqOp<point>(),
            -greatPoint
        );
        syncTools::syncPointPositions
        (
            fvMesh_,
            p2,
            minEqOp<point>(),
            greatPoint
        );
        cPts = (p1 + p2)/2.0;
    }

    forAll(fixedVelocity_, pointI)
    {
        if (fixedVelocity_[pointI] != vector(GREAT, GREAT, GREAT))
        {
            cPts[pointI] =
                fvMesh_.points()[pointI]
              + fvMesh_.time().deltaTValue()*fixedVelocity_[pointI];
        }
    }

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::velocityLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();

    fvVectorMatrix UEqn
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellMotionU_,
            "laplacian(diffusivity,cellMotionU)"
        )
    );

    fixedVelocity_ = vector(GREAT, GREAT, GREAT);

    boolList markedCells(fvMesh_.nCells(), false);
    DynamicList<label> cells(fvMesh_.nCells());
    DynamicList<vector> fixedValues(fvMesh_.nCells());

    forAll(pointMotionU_.boundaryField(), patchI)
    {
        labelList tempCells(0);
        vectorField tempFixedValues(0);

        // In most of boundaries this isn't implemented hance will do nothing
        pointMotionU_.boundaryFieldRef()[patchI].manipulateMatrix
        (
            tempCells,
            tempFixedValues
        );
        forAll(tempCells, i)
        {
            const label celli = tempCells[i];

            if (!markedCells[celli])
            {
                markedCells[celli] = true;
                cells.append(celli);
                fixedValues.append(tempFixedValues[i]);
            }
        }

        // In most of boundaries this isn't implemented hance will do nothing
        pointMotionU_.boundaryFieldRef()[patchI].setField(fixedVelocity_);
    }

    cells.shrink();
    fixedValues.shrink();

    if (cells.size())
    {
        UEqn.setValues(cells, vectorField(fixedValues));
    }

    UEqn.solveSegregatedOrCoupled(UEqn.solverDict());
}


//void Foam::velocityLaplacianFvMotionSolver::movePoints(const pointField& p)
//{
//    // Movement of pointMesh and volPointInterpolation already
//    // done by polyMesh,fvMesh
//}


void Foam::velocityLaplacianFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    velocityMotionSolver::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


void Foam::velocityLaplacianFvMotionSolver::mapMesh(const polyMeshMap& map)
{
    velocityMotionSolver::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
