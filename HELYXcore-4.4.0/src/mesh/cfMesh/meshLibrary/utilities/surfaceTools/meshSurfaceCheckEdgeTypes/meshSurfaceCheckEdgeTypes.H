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
    meshSurfaceCheckEdgeTypes

Description
    Checks if any of the surface is tangled at any of its vertices. This
    is performed by checking the dot product between the vertex normal and
    the normals of the faces attached to the vertex.

SourceFiles
    meshSurfaceCheckEdgeTypes.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCheckEdgeTypes_H
#define meshSurfaceCheckEdgeTypes_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/Time/Time.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
            Class meshSurfaceCheckEdgeTypes Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceCheckEdgeTypes
{
    // Private data

        //- Mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- A list which classifies each edge
        List<direction> edgeType_;


    // Private member functions

        //- Check feature edges and classify them
        void classifyEdges();


public:

    // Enumerators

        enum edgeTypes_
        {
            NONE = 0,
            PATCHEDGE = 1,
            FEATUREEDGE = 2,
            CONVEXEDGE = 4,
            CONCAVEEDGE = 8,
            UNDETERMINED = 16
        };

    // Constructors

        //- Construct from meshSurfaceEngine
        meshSurfaceCheckEdgeTypes(const meshSurfaceEngine&);

        //- Disallow default bitwise copy construction
        meshSurfaceCheckEdgeTypes(const meshSurfaceCheckEdgeTypes&) = delete;


    //- Destructor
    ~meshSurfaceCheckEdgeTypes();


    // Member Functions

        //- Return the information about classified edges
        const List<direction>& edgeTypes() const
        {
            return edgeType_;
        }

        //- Return indices of convex edges
        void convexEdges(labelLongList& convexEdges) const;

        //- Return indices of concave edges
        void concaveEdges(labelLongList& concaveEdges) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceCheckEdgeTypes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
