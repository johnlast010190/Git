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
    tetMeshExtractorOctree

Description
    A class which extracts tet mesh out of an octree structure

SourceFiles
    tetMeshExtractorOctree.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshExtractorOctree_H
#define tetMeshExtractorOctree_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "utilities/tetrahedra/tetCreatorOctree/tetCreatorOctree.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tetMeshExtractorOctree Declaration
\*---------------------------------------------------------------------------*/

class tetMeshExtractorOctree
{
    // Private data

        //- Create tets
        tetCreatorOctree tetCreator_;

        //- Reference to the mesh
        polyMeshGen& mesh_;


    // Private member functions

        //- Copy tetPoints_ into mesh
        void createPoints();

        //- Create mesh data
        void createPolyMesh();


public:

    // Constructors

        //- Construct from octree and mesh data
        tetMeshExtractorOctree
        (
            const meshOctree& octree,
            const IOdictionary& meshDict,
            polyMeshGen& mesh
        );

        //- Disallow default bitwise copy construction
        tetMeshExtractorOctree (const tetMeshExtractorOctree&) = delete;


    //- Destructor
    ~tetMeshExtractorOctree();


    // Member Functions

        //- Create mesh
        void createMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tetMeshExtractorOctree&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
