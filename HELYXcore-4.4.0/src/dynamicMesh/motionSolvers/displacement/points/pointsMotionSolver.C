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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/points/pointsMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointsMotionSolver, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::pointVectorField Foam::pointsMotionSolver::readPoints0
(
    const polyMesh& mesh
)
{
    const word instance =
        mesh.time().findInstance
        (
            mesh.meshDir(),
            "points0",
            IOobject::READ_IF_PRESENT
        );

    if (instance != mesh.time().constant())
    {
        // points0 written to a time folder

        // backward compatibility check for v4.2 and earlier
        // first read IOobject
        IOobject ioPoints0
        (
            "points0",
            instance,
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        // then check its type
        if (ioPoints0.typeHeaderOk<pointIOField>(true, true, false))
        {
            DeprecationWarningInFunction
            (
                "pointIOField",
                "field type for points0",
                40300,
                "Please use pointVectorField instead."
            );

            Info<< "    Converting points0 from pointIOField"
                 << " to pointVectorField." << endl;
            pointIOField points0IO(ioPoints0);

            ioPoints0.readOpt() = IOobject::NO_READ;
            pointVectorField points0
            (
                ioPoints0,
                pointMesh::New(mesh),
                dimensionedVector(dimLength, Zero)
            );
            points0.primitiveFieldRef() = points0IO;

            return points0;
        }

        return pointVectorField
        (
            ioPoints0,
            pointMesh::New(mesh)
        );
    }
    else
    {
        // Return copy of original mesh points

        pointIOField points
        (
            IOobject
            (
                "points",
                mesh.time().findInstance(mesh.meshDir(), "points"),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        pointVectorField points0
        (
            IOobject
            (
                "points",
                instance,
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh),
            dimensionedVector("points0", dimLength, Zero)
        );

        points0.primitiveFieldRef() = points;

        return points0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointsMotionSolver::pointsMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(readPoints0(mesh))
{
    if
    (
        dict.lookupOrDefault<Switch>("resetPoints0", false)
     || (!dict.found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "resetting points0" << endl;
        points0_.primitiveFieldRef() = mesh.points();
    }

    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file "
            <<  typeFilePath<pointIOField>
                (
                    IOobject
                    (
                        "points",
                        mesh.time().constant(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            << exit(FatalError);
    }

    // If points0 was obtained from latest mesh points, write it out to
    // ensure stable restart behaviour.
    if (points0_.name() == "points")
    {
        fileName pointsFileName(points0_.objectPath());
        points0_.rename("points0");
        // Don't write into constant to avoid confusion
        points0_.instance() = mesh.time().timeName();
        points0_.write();
    }
}


Foam::pointsMotionSolver::pointsMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const pointVectorField& points0,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    points0_(points0)
{
    if
    (
        dict.lookupOrDefault<Switch>("resetPoints0", false)
     || (!dict.found("resetPoints0") && isIncrementalMotion())
    )
    {
        Info<< "resetting points0" << endl;
        points0_.primitiveFieldRef() = mesh.points();
    }

    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file " << points0.filePath()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointsMotionSolver::~pointsMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointsMotionSolver::movePoints(const pointField&)
{}


void Foam::pointsMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    NotImplemented;
}


void Foam::pointsMotionSolver::mapMesh(const polyMeshMap& map)
{
    points0_.primitiveFieldRef() = mesh().points();
}


void Foam::pointsMotionSolver::distribute
(
    const polyDistributionMap& map
)
{}


bool Foam::pointsMotionSolver::write() const
{
    if (dynamic_cast<const fvMesh&>(mesh()).topoChanging())
    {
        points0_.instance() = mesh().time().timeName();
        points0_.write();
    }

    return motionSolver::write();
}


// ************************************************************************* //
