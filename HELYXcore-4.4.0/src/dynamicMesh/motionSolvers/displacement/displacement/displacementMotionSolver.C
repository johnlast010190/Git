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
    (c) 2016 OpenCFD Ltd.
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMotionSolver, 0);
    defineRunTimeSelectionTable(displacementMotionSolver, displacement);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMotionSolver::displacementMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    pointsMotionSolver(mesh, dict, type),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


Foam::displacementMotionSolver::displacementMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointVectorField& points0,
    const word& type
)
:
    pointsMotionSolver(mesh, dict, points0, type),
    pointDisplacement_
    (
        IOobject(pointDisplacement, "pointDisplacement"),
        pointDisplacement
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::displacementMotionSolver>
Foam::displacementMotionSolver::New
(
    const word& solverTypeName,
    const polyMesh& mesh,
    const dictionary& solverDict,
    const pointVectorField& pointDisplacement,
    const pointVectorField& points0
)
{
    Info<< "Selecting motion solver: " << solverTypeName << endl;

    libs.open(solverDict, "motionSolverLibs");

    const auto ctor =
        ctorTableLookup
        (
            "solver type",
            displacementConstructorTable_(),
            solverTypeName
        );
    return autoPtr<displacementMotionSolver>
    (
        ctor
        (
            mesh,
            solverDict,
            pointDisplacement,
            points0
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementMotionSolver::~displacementMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementMotionSolver::mapMesh(const polyMeshMap& map)
{
    pointsMotionSolver::mapMesh(map);
    pointDisplacement_.forceAssign(vector::zero);
}


// ************************************************************************* //
