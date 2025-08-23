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
    cartesianMeshExtractor

Description
    Generates the cartesian mesh template from the octree

SourceFiles
    cartesianMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshExtractor_H
#define cartesianMeshExtractor_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOdictionary;

/*---------------------------------------------------------------------------*\
                    Class cartesianMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshExtractor
{
    // Private data

        //- Reference to the octree addressing
        meshOctreeAddressing octreeCheck_;

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- Decompose split hex cells
        bool decomposeSplitHexes_;

        //- Cell label for a given leaf
        labelList* leafCellLabelPtr_;


    // Private member functions

        //- Delete all freestore data
        void clearOut();

        //- Create vertices and pointLeaves addressing
        void createPointsAndAddressing();

        //- Create mesh data
        void createPolyMesh();

        //- Decompose split hexes into pyramids and tets
        void decomposeSplitHexesIntoTetsAndPyramids();


public:

    // Constructors

        //- Construct from octree and mesh data
        cartesianMeshExtractor
        (
            meshOctree& octree,
            const IOdictionary& meshDict,
            polyMeshGen& mesh
        );

        //- Disallow default bitwise copy construction
        cartesianMeshExtractor(const cartesianMeshExtractor&) = delete;


    // Destructor
    ~cartesianMeshExtractor();


    // Member Functions

        //- Decompose split hexes into standard cells
        void decomposeSplitHexes();

        //- Create the mesh with the above options
        void createMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cartesianMeshExtractor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
