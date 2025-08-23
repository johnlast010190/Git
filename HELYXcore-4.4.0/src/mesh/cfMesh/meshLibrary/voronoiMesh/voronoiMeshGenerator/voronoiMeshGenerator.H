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
    voronoiMeshGenerator

Description
    Creates voronoi mesh template from the octree

SourceFiles
    voronoiMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef voronoiMeshGenerator_H
#define voronoiMeshGenerator_H

#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "utilities/workflowControls/workflowControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;
class meshOctree;
class Time;

/*---------------------------------------------------------------------------*\
                     Class voronoiMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class voronoiMeshGenerator
{
    // Private data

        //- Reference to Time
        const Time& runTime_;

        //- Pointer to the surface
        const triSurf* surfacePtr_;

        //- Pointer to the modifier surface
        const triSurf* modSurfacePtr_;

        //- Pointer to the octree
        meshOctree* octreePtr_;

        //- Pointer to the list patches for boundary vertices
        labelList* pointRegionsPtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- Mesh
        polyMeshGen mesh_;

        //- Workflow controller
        workflowControls controller_;


    // Private member functions

        //- Create voronoi mesh
        void createVoronoiMesh();

        //- Prepare mesh surface
        void surfacePreparation();

        //- Map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- Assign patches to boundary faces
        void extractPatches();

        //- Capture edges and corners
        void mapEdgesAndCorners();

        //- Optimise surface mesh
        void optimiseMeshSurface();

        //- Add boundary layers
        void generateBoudaryLayers();

        //- Mesh optimisation
        void optimiseFinalMesh();

        //- Re-project points back on the surface mesh after back-scaling
        void projectSurfaceAfterBackScaling();

        //- Refine boundary layer
        void refBoundaryLayers();

        //- Replace boundaries
        void replaceBoundaries();

        //- Renumber the mesh
        void renumberMesh();

        //- Generate mesh
        void generateMesh();


public:

    // Constructors

        //- Construct from time
        voronoiMeshGenerator(const Time&);

        //- Disallow default bitwise copy construction
        voronoiMeshGenerator(const voronoiMeshGenerator&) = delete;


    //- Destructor
    ~voronoiMeshGenerator();


    // Member Functions

        //- write the mesh
        void writeMesh() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const voronoiMeshGenerator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
