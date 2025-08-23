
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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"
#include "patchToPatch/patchToPatch/patchToPatch.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("source");
    argList::validArgs.append("target");
    argList::validArgs.append("method");

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createPolyMesh.H"

    const polyPatch& srcPatch = mesh.boundaryMesh()[args[1]];
    const polyPatch& tgtPatch = mesh.boundaryMesh()[args[2]];
    const word& method = args[3];

    cpuTime time;

    /*
    AMIInterpolation
    (
        refCast<const directPolyPatch>(srcPatch),
        refCast<const directPolyPatch>(tgtPatch),
        faceAreaIntersect::tmMesh,
        true,
        "faceAreaWeightAMI"
    );

    Info<< nl << "AMI" << ": Completed in "
        << time.cpuTimeIncrement() << " s" << nl << endl;
    */

    patchToPatch::New(method, false)->update
    (
        refCast<const directPolyPatch>(srcPatch),
        srcPatch.pointNormals(),
        refCast<const directPolyPatch>(tgtPatch)
    );

    Info<< nl << patchToPatch::typeName << ": Completed in "
        << time.cpuTimeIncrement() << " s" << nl << endl;

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
