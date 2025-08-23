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

\*---------------------------------------------------------------------------*/

#include "include/demandDrivenData.H"
#include "utilities/smoothers/geometry/meshOptimizer/meshOptimizer.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizer.H"
#include "utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapper.H"
#include "utilities/octrees/meshOctree/meshOctree.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::optimizeSurface(const meshOctree& octree)
{
    Info<< "Optimizing positions of surface nodes" << endl;

    meshSurfaceEngine& mse = const_cast<meshSurfaceEngine&>(meshSurface());
    meshSurfaceOptimizer surfaceOptimizer(mse, octree);

    if (enforceConstraints_)
        surfaceOptimizer.enforceConstraints(badPointsSubsetName_);

    surfaceOptimizer.optimizeSurface();

    meshSurfaceMapper(mse, octree).mapVerticesOntoSurfacePatches();

    clearSurface();

    Info<< "Finished optimizing positions of surface nodes" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
