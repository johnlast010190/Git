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
    (c) 2013 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMotionSolvers/displacement/laplacianPeriodic/displacementLaplacianPeriodicFvMotionSolver.H"
#include "motionInterpolation/motionInterpolation/motionInterpolation.H"
#include "motionDiffusivity/motionDiffusivity/motionDiffusivity.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshTools/meshTools.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/pointPatchFields/basic/value/valuePointPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLaplacianPeriodicFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLaplacianPeriodicFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        displacementLaplacianPeriodicFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::displacementLaplacianPeriodicFvMotionSolver::
createMaxDisplacementField()
{
    IOobject maxPDispHeader
    (
        "maxPointDisplacement",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (maxPDispHeader.headerOk())
    {
        startFromDisp_ = true;
    }

    // Read maxPointDisplacement if present
    maxPointDisplacement_.reset
    (
        new pointVectorField
        (
            maxPDispHeader,
            pointDisplacement_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLaplacianPeriodicFvMotionSolver::
displacementLaplacianPeriodicFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementLaplacianFvMotionSolver(mesh, dict),
    amplification_(Function1<scalar>::New("amplification", coeffDict())),
    maxDisplacementTime_
    (
        readScalar(coeffDict().lookup("maxDisplacementTime"))
    ),
    startFromDisp_(false)
{
    createMaxDisplacementField();
}


Foam::displacementLaplacianPeriodicFvMotionSolver::
displacementLaplacianPeriodicFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointVectorField& points0
)
:
    displacementLaplacianFvMotionSolver
    (
        mesh, dict, pointDisplacement, points0
    ),
    amplification_(Function1<scalar>::New("amplification", coeffDict())),
    maxDisplacementTime_
    (
        readScalar(coeffDict().lookup("maxDisplacementTime"))
    ),
    startFromDisp_(false)
{
    createMaxDisplacementField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLaplacianPeriodicFvMotionSolver::
~displacementLaplacianPeriodicFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementLaplacianPeriodicFvMotionSolver::
computeCellDisplacement()
{
    displacementLaplacianFvMotionSolver::solve();
}


void Foam::displacementLaplacianPeriodicFvMotionSolver::
computePointDisplacement()
{
    // Compute cell displacement
    computeCellDisplacement();

    // Do vol-to-point interpolation
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    // Initialise current point field
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    pointField& curPoints = tcurPoints.ref();

    // Constrain the given field (2D-correction + fix buffer)
    constrainInternalPoints(curPoints);

    // Update pointDisplacement after correction
    pointDisplacement_.primitiveFieldRef() = curPoints-points0();

    // Store point displacement (=max since called at maxDisplacementTime_)
    maxPointDisplacement_->primitiveFieldRef() = curPoints-points0();
    forAll(pointDisplacement_.boundaryField(), pI)
    {
        maxPointDisplacement_->boundaryFieldRef()[pI].forceAssign
        (
            pointDisplacement_.boundaryField()[pI]
        );
    }
}


Foam::tmp<Foam::pointField>
Foam::displacementLaplacianPeriodicFvMotionSolver::curPoints() const
{
    bool finalCorrector =
        fvMesh_.time().lookupObject<helyxSolve>(helyxSolve::typeName)
        .finalCorrector("initCorrector");

    // At finalCorrector set mesh to actual position
    if (fvMesh_.time().timeIndex() == 0 && !finalCorrector)
    {
        return curPoints(maxDisplacementTime_);
    }

    return curPoints(fvMesh_.time().value());
}


Foam::tmp<Foam::pointField>
Foam::displacementLaplacianPeriodicFvMotionSolver::curPoints
(
    const scalar& time
) const
{
    // Compute internal current points based on max displacement and
    // amplification
    tmp<pointField> tcurPoints
    (
        points0()
      + maxPointDisplacement_->primitiveField()
      * amplification_->value(time)
    );

    return tcurPoints;
}


void Foam::displacementLaplacianPeriodicFvMotionSolver::solve()
{
    // At first call compute max displacement (boundary and internal)
    // else do nothing
    if (fvMesh_.time().timeIndex() == 0 && !startFromDisp_)
    {
        // Set time to maxDisplacementTime_
        Time& runTime = const_cast<Time&>(fvMesh_.time());

        scalar currentTime = runTime.value();
        scalar deltaT = runTime.deltaTValue();

        runTime.setDeltaT(maxDisplacementTime_-currentTime);
        runTime.setTime(maxDisplacementTime_, runTime.timeIndex());

        // Compute maxPointDisplacement_ field inclusive boundaries
        computePointDisplacement();

        // Re-set time to current
        runTime.setDeltaT(deltaT);
        runTime.setTime(currentTime, runTime.timeIndex());
    }
}


void Foam::displacementLaplacianPeriodicFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    NotImplemented;
}


void Foam::displacementLaplacianPeriodicFvMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    NotImplemented;
}


// ************************************************************************* //
