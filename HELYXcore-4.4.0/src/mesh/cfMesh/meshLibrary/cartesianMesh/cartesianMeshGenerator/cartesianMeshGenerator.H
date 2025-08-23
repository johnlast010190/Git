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
    cartesianMeshGenerator

Description
    Creates cartesian mesh from the octree

SourceFiles
    cartesianMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshGenerator_H
#define cartesianMeshGenerator_H

#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "utilities/workflowControls/workflowControls.H"
//#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;
class meshOctree;
class Time;

/*---------------------------------------------------------------------------*\
                           Class cartesianMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshGenerator
{
    // Private data

        //- Reference to Time
        const Time& db_;

        //- Pointer to the surface
        const triSurf* surfacePtr_;

        //- Pointer to the modified surface
        const triSurf* modSurfacePtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- Pointer to the octree
        meshOctree* octreePtr_;

        //- Mesh
        polyMeshGen mesh_;

        //- Workflow controller
        workflowControls controller_;


    // Private member functions

        //- Create cartesian mesh
        void createCartesianMesh();

        //- Prepare mesh surface
        void surfacePreparation();

        //- Map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- Capture edges and corners
        void extractPatches();

        //- Map points at feature edges and corners
        void mapEdgesAndCorners();

        //- Optimise surface mesh
        void optimiseMeshSurface();

        //- Add boundary layers
        void generateBoundaryLayers();

        //- Refine boundary layers
        void refBoundaryLayers();

        //- Mesh optimisation
        void optimiseFinalMesh();

        //- Re-project points back on the surface mesh after back-scaling
        void projectSurfaceAfterBackScaling();

        //- Replace boundaries
        void replaceBoundaries();

        //- Renumber the mesh
        void renumberMesh();

        //- Generate mesh
        void generateMesh();


public:

    // Constructors

        //- Construct from time
        cartesianMeshGenerator(const Time&);

        //- Disallow default bitwise copy construction
        cartesianMeshGenerator(const cartesianMeshGenerator&) = delete;


    //- Destructor
    ~cartesianMeshGenerator();


    // Member Functions

        //- Write the mesh
        void writeMesh() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cartesianMeshGenerator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
