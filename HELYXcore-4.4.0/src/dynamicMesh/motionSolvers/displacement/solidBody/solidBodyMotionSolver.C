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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/solidBody/solidBodyMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/Fields/transformField/transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        solidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionSolver::solidBodyMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointsMotionSolver(mesh, dict, typeName),
    transform_(transformation())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionSolver::~solidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::solidBodyMotionSolver::curPoints() const
{
    transform_ = transformation();

    if (moveAllPoints())
    {
        return transformPoints(transform_, points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, setPointIndices()) =
            transformPoints
            (
                transform_,
                pointField(points0_, setPointIndices())
            );

        return ttransformedPts;
    }
}


void Foam::solidBodyMotionSolver::solve()
{
    if (mesh().changing())
    {
        Info<< "Update zones" << endl;
        updateSetPointIndices();
    }

    if
    (
        coeffDict().lookupOrDefault<Switch>("resetPoints0", false)
     || (!coeffDict().found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "Resetting points0" << endl;
        points0() = mesh().points();
    }
}


void Foam::solidBodyMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    // pointMesh already updates pointFields
    //motionSolver::topoChange(map);

    // Map points0_. Bit special since we somehow have to come up with
    // a sensible points0 position for introduced points.
    // Find out scaling between points0 and current points

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        map.hasMotionPoints()
      ? map.preMotionPoints()
      : mesh().points()
    );

    pointField newPoints0(map.pointMap().size());

    forAll(newPoints0, pointi)
    {
        label oldPointi = map.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            label masterPointi = map.reversePointMap()[oldPointi];

            if (masterPointi == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
            else
            {
                newPoints0[pointi] =
                    transform_.invTransformPoint(points[pointi]);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot determine co-ordinates of introduced vertices."
                << " New vertex " << pointi << " at co-ordinate "
                << points[pointi] << exit(FatalError);
        }
    }

    twoDCorrectPoints(newPoints0);

    points0_.transfer(newPoints0);
}


// ************************************************************************* //
