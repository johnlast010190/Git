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
    meshSurfaceCheckInvertedVertices

Description
    Checks if any of the surface is tangled at any of its vertices. This
    is performed by checking the dot product between the vertex normal and
    the normals of the faces attached to the vertex.

SourceFiles
    meshSurfaceCheckInvertedVertices.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCheckInvertedVertices_H
#define meshSurfaceCheckInvertedVertices_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/Time/Time.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfacePartitioner;

/*---------------------------------------------------------------------------*\
            Class meshSurfaceCheckInvertedVertices Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceCheckInvertedVertices
{
    // Private data

        //- Mesh surface partitioner
        const meshSurfacePartitioner& surfacePartitioner_;

        //- Active surface points
        const boolList* activePointsPtr_;

        //- Set of inverted vertices
        labelHashSet invertedVertices_;


    // Private member functions

        //- Check vertices by calculating dot products
        void checkVertices();


public:

    // Constructors

        //- Construct from meshSurfacePartitioner
        meshSurfaceCheckInvertedVertices
        (
            const meshSurfacePartitioner&
        );

        //- Construct from meshSurfacePartitioner and a list of active points
        meshSurfaceCheckInvertedVertices
        (
            const meshSurfacePartitioner&,
            const boolList& activePoints
        );

        //- Disallow default bitwise copy construction
        meshSurfaceCheckInvertedVertices
        (
            const meshSurfaceCheckInvertedVertices&
        ) = delete;


    //- Destructor
    ~meshSurfaceCheckInvertedVertices();


    // Member Functions

        //- Return the labels of inverted vertices
        inline const labelHashSet& invertedVertices() const
        {
            return invertedVertices_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceCheckInvertedVertices&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
