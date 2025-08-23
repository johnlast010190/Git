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
    (c) 2011-2016 OpenFOAM Foundation

Description
    Calculate and write the distance-to-wall field for a moving mesh.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;

    Info<< "Time now = " << runTime.timeName() << endl;

    // Wall-reflection vectors
    const volVectorField& n = wallDist::New(mesh).n();
    n.write();

    // Wall distance
    const volScalarField& y = wallDist::New(mesh).y();
    y.write();

    runTime++;

    Info<< "Time now = " << runTime.timeName() << endl;

    // Move points

    boundBox meshBb(mesh.points());

    pointField newPoints(mesh.points());

    const point half = meshBb.midpoint();

    forAll(newPoints, pointi)
    {
        point& pt = newPoints[pointi];

        // expand around half
        pt.y() += pt.y() - half.y();
    }

    mesh.movePoints(newPoints);
    mesh.write();
    n.write();
    y.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
