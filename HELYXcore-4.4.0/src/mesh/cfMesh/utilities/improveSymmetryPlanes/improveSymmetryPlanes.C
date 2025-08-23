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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description
    Ensures that all mesh points belonging to a symmetryPlane are
    in a plane.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/smoothers/geometry/meshOptimizer/symmetryPlaneOptimisation/symmetryPlaneOptimisation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "include/setRootCase.H"
#include "include/createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    Info<< "Starting optimisation of symmetry planes" << endl;
    symmetryPlaneOptimisation(pmg).optimizeSymmetryPlanes();

    Info<< "Writing mesh" << endl;
    pmg.write();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
