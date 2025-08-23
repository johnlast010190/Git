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
    (c) 2020-2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    meshOptimizer

Description
    Mesh smoothing without any topological changes

SourceFiles
    meshOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef meshOptimizer_H
#define meshOptimizer_H

#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "primitives/bools/lists/boolList.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "utilities/containers/DynList/DynList.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "meshes/boundBox/boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class plane;
class partTetMesh;

/*---------------------------------------------------------------------------*\
                        Class meshOptimizer Declaration
\*---------------------------------------------------------------------------*/

class meshOptimizer
{
    // Private data

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- Location of vertex (internal, boundary, edge, corner)
        List<direction> vertexLocation_;

        //- Locked faces which shall not be changed
        labelLongList lockedFaces_;

        //- Mesh surface
        mutable meshSurfaceEngine* msePtr_;

        //- Enforce constraints
        bool enforceConstraints_;

        //- Name of the subset contaning tangled points
        word badPointsSubsetName_;


    // Private member functions

        //- Return mesh surface
        const meshSurfaceEngine& meshSurface() const;
        void clearSurface();

        //- Find problematic faces
        label findBadFaces(labelHashSet&, const boolList&) const;
        label findLowQualityFaces(labelHashSet&, const boolList&) const;

        //- Mark point locations
        void calculatePointLocations();

    // Nested classes

        class laplaceSmoother
        {
            // Private data

                //- Reference to the mesh
                polyMeshGen& mesh_;

                //- Location of vertex (internal, boundary, edge, corner)
                const List<direction>& vertexLocation_;


            // Private member functions

                //- Smooth the node using the laplacian smoother
                //  new position is the average of the neighbouring vertices
                void laplacian(const labelLongList&, const label);
                void laplacianSurface(const labelLongList&, const label);

                void laplacianParallel
                (
                    const labelLongList& procPoints,
                    const bool smoothOnlySurfaceNodes = false
                );

                //- Smooth the node using the laplacian smoother
                //  new position is the average of the centres of faces attached
                //  to the vertex
                void laplacianPC(const labelLongList&, const label);
                void laplacianPCParallel(const labelLongList& procPoints);

                //- Smooth the node using the laplacian smoother
                //  new position is the average of the centres of faces attached
                //  to the vertex
                void laplacianWPC(const labelLongList&, const label);
                void laplacianWPCParallel(const labelLongList& procPoints);

                //- Update geometry after smoothing
                void topoChangeGeometry(const labelLongList& smoothPoints);


        public:

            // Constructor

                //- Construct from mesh and vertex locations
                laplaceSmoother(polyMeshGen&, const List<direction>&);

                //- Disallow default bitwise copy construction
                laplaceSmoother(const laplaceSmoother&) = delete;


            //- Destructor
            ~laplaceSmoother();


            // Member Functions

                //- New position is the average of the neighbouring vertices
                void optimizeLaplacian(const label nIterations = 1);
                void optimizeLaplacian
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );

                //- New position of surface point is the average of
                //  the neighbouring surface vertices
                void optimizeSurfaceLaplacian
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );

                //- New positions are the average of the centres of the cells
                //  adjacent to the vertex
                void optimizeLaplacianPC(const label nIterations = 1);
                void optimizeLaplacianPC
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );

                //- New positions are the average of the centres of the cells
                //  adjacent to the vertex weighted by cell volumes
                void optimizeLaplacianWPC(const label nIterations = 1);
                void optimizeLaplacianWPC
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );


            // Member Operators

                //- Disallow default bitwise assignment
                void operator=(const laplaceSmoother&) = delete;
        };


    // Enumerators
        enum vertexType_
        {
            INSIDE = 1,
            BOUNDARY = 2,
            EDGE = 4,
            CORNER = 8,
            PARALLELBOUNDARY = 16,
            LOCKED = 32
        };

public:

    // Constructors

        //- Construct from mesh
        meshOptimizer(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        meshOptimizer(const meshOptimizer&) = delete;


    //- Destructor
    ~meshOptimizer();


    // Member Functions

        //- Flag stopping the meshing process if it is not possible
        //  to untangle the surface without sacrificing geometry constraints
        //  Points which must be moved away from the required position are
        //  stored into a point subset
        void enforceConstraints(const word subsetName="badPoints");

        //- Lock cells which are error free
        void lockErrorFreeCells
        (
            const label nBufferLayers,
            const scalar minFaceArea
        );

        //- Lock the cells which shall not be modified
        template<class labelListType>
        inline void lockCells(const labelListType&);

        //- Lock cells which shall not be modified
        //  contained inside the specified subset
        void lockCellsInSubset(const word& subsetName);

        //- Lock the faces whih shall not be modified
        template<class labelListType>
        inline void lockFaces(const labelListType&);

        //- Lock faces which shall not be modified
        //  stored in a face subset
        void lockFacesInSubset(const word& subsetName);

        //- Lock points which shall not be moved
        template<class labelListType>
        inline void lockPoints(const labelListType&);

        //- Lock points which shall not be modified
        //  given as a point subset
        void lockPointsInSubset(const word& subsetName);

        //- Reset to default constraints
        void removeUserConstraints();

        //- Smooth surface vertices
        void optimizeSurface(const meshOctree&);

        //- Performs mesh untangling based on detected negative normals
        //  a global iteration consists of maxNumIterations iterations when
        //  only internal points are moved and maxNumSurfaceIterations when
        //  surface points are moved
        void untangleMeshFV
        (
            const label maxNumGlobalIterations = 10,
            const label maxNumIterations = 50,
            const label maxNumSurfaceIterations = 2,
            const label errorBufferLayers = -1,
            const bool relaxedCheck = false,
            const bool relaxedBoundaryCheck = false,
            const bool checkWarped = true,
            const scalar minFaceArea = VSMALL
        );

        //- Performs optimisation of boundary layer cells
        //  it is applied to the first boundary layer, only
        void optimizeBoundaryLayer(const bool addBufferLayer = true);

        //- Performs untangling of boundary layer cells
        //  if the user requests for it
        void untangleBoundaryLayer();

        //- Performs mesh optimisation for faces with high non-orthogonality
        //  and skewness
        void optimizeLowQualityFaces(const label maxNumIterations = 10);

        //- Perform optimization of cells near the boundary
        //  the use shall specify the number of iterations and the number
        //  of cell layers which shall be taken into account
        void optimizeMeshNearBoundaries
        (
            const label maxNumIterations = 10,
            const label numLayersOfCells = 2
        );

        //- Final optimisation
        void optimizeMeshFV
        (
            const label numLaplaceIterations = 5,
            const label maxNumGlobalIterations = 10,
            const label maxNumIterations = 50,
            const label maxNumSurfaceIterations = 2,
            const label errorBufferLayers = -1,
            const bool relaxedCheck = false,
            const bool relaxedBoundaryCheck = false,
            const bool checkWarped = true,
            const scalar minFaceArea = VSMALL
        );

        //- Greedy optimisation until the mesh can be improved
        void optimizeMeshFVBestQuality
        (
            const label maxNumIterations = 50,
            const scalar threshold = 0.1
        );

        //- See if mesh orthogonlaity can be improved
        void optimizeMeshFVOrthogonality
        (
            const label maxNumIterations = 50,
            const scalar threshold = 0.342
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshOptimizer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/smoothers/geometry/meshOptimizer/meshOptimizerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
