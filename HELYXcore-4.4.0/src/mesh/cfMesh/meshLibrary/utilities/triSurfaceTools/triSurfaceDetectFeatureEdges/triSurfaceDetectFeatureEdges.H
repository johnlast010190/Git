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
    triSurfaceDetectFeatureEdges

Description
    Divides the surface mesh into regions bounded by feature edges

SourceFiles
    triSurfaceDetectFeatureEdges.C
    triSurfaceDetectFeatureEdgesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDetectFeatureEdges_H
#define triSurfaceDetectFeatureEdges_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceDetectFeatureEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDetectFeatureEdges
{
    // Private data

        //- Reference to triSurf
        triSurf& surf_;

        //- Detected feature edges
        List<direction> featureEdges_;

        //- Angle tolerance
        scalar angleTolerance_;


    // Private member functions

        //- Detect feature edges by check angle between the normals
        //  of the two neighbouring facets
        void detectFeatureEdgesAngleCriterion();


public:

    // Constructors

        //- Construct from triSurface and tolerance
        triSurfaceDetectFeatureEdges
        (
            triSurf& surface,
            const scalar angleDeviation = 45.0
        );

        //- Disallow default bitwise copy construction
        triSurfaceDetectFeatureEdges
        (
            const triSurfaceDetectFeatureEdges&
        ) = delete;


    //- Destructor
    ~triSurfaceDetectFeatureEdges();


    // Member Functions

        //- Store regions in a graph
        void detectFeatureEdges();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceDetectFeatureEdges&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
