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
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    triSurfaceExtrude2DEdges

Description
    Extrudes edges in x-y plane into a triangulation used for 2D meshing

SourceFiles
    triSurfaceExtrude2DEdges.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceExtrude2DEdges_H
#define triSurfaceExtrude2DEdges_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceExtrude2DEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceExtrude2DEdges
{
    // Private data

        //- reference to triSurf
        const triSurf& surf_;


public:

    // Constructors

        //- Construct from octree
        triSurfaceExtrude2DEdges(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfaceExtrude2DEdges(const triSurfaceExtrude2DEdges&) = delete;


    // Destructor
    ~triSurfaceExtrude2DEdges();


    // Member Functions

        //- extrudes edges into a triangulation
        void extrudeSurface(triSurf&) const;

        const triSurf* extrudeSurface() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceExtrude2DEdges&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
