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
    (c) 2023 FOSS GP
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/modes/modesMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(modesMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        modesMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::modesMotionSolver::modesMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    fileNames_(coeffDict().lookup<wordReList>("files")),
    amplifications_(fileNames_.size()),
    modes_(fileNames_.size()),
    points0_(mesh.points()),
    appendToCurrentPoints_
        (coeffDict().lookupOrDefault<bool>("appendToCurrentPoints", true))
{
    forAll(fileNames_, fi)
    {
        modes_.set
        (
            fi,
            new pointVectorField
            (
                IOobject
                (
                    fileNames_[fi],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pointMesh::New(mesh)
            )
        );
    }

    forAll(fileNames_, fi)
    {
        amplifications_.set
        (
            fi, Function1<scalar>::New(fileNames_[fi], coeffDict())
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::modesMotionSolver::curPoints() const
{
    auto tcurPoints
    (
        tmp<pointField>::New
        (
            appendToCurrentPoints_
          ? mesh().points()
          : points0_
        )
    );
    pointField& curPoints = tcurPoints.ref();
    forAll(modes_, mi)
    {
        curPoints +=
            modes_[mi]*amplifications_[mi].value(mesh().time().value());
    }

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::modesMotionSolver::movePoints(const pointField& p)
{
    // No local data that needs adapting.
}


void Foam::modesMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    // pointMesh already updates pointFields.

    motionSolver::topoChange(map);
}


void Foam::modesMotionSolver::mapMesh(const polyMeshMap& map)
{
    // pointMesh already updates pointFields.
    //motionSolver::mapMesh(map);
}


// ************************************************************************* //
