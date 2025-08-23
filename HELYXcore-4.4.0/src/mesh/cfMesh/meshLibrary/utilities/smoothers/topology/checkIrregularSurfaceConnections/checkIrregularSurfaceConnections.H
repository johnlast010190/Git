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
    checkIrregularSurfaceConnections

Description
    Check mesh surface for for connections which cannot exist in a valid mesh.
    In case some irregular connections exist, it is possible to fix them.

SourceFiles
    checkIrregularSurfaceConnections.C

\*---------------------------------------------------------------------------*/

#ifndef checkIrregularSurfaceConnections_H
#define checkIrregularSurfaceConnections_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class checkIrregularSurfaceConnections Declaration
\*---------------------------------------------------------------------------*/

class checkIrregularSurfaceConnections
{
    // Private data

        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- Pointer to meshSurfaceEngine
        mutable meshSurfaceEngine* meshSurfacePtr_;


    // Private member functions

        //- Create meshSurfaceEngine
        inline const meshSurfaceEngine& surfaceEngine() const
        {
            if (!meshSurfacePtr_)
                meshSurfacePtr_ = new meshSurfaceEngine(mesh_);

            return *meshSurfacePtr_;
        }

        //- Clear meshSurfaceEngine. Used in case if the mesh is modified
        inline void clearMeshEngine()
        {
            deleteDemandDrivenData(meshSurfacePtr_);
        }

        //- Check if there exist vertices with more than one groups of cells
        //  attached to it
        bool checkAndFixCellGroupsAtBndVertices
        (
            labelHashSet& badVertices,
            const bool removeConnections = false
        );

        //- Check if all boundary edges are connected
        //  to two boundary faces, only!
        bool checkEdgeFaceConnections
        (
            labelHashSet& badVertices,
            const bool removeCells = false
        );

        //- Check if all boundary face connected to a boundary vertex
        //  are connected into a single group of faces that can be visited
        //  over edges
        bool checkFaceGroupsAtBndVertices
        (
            labelHashSet& badVertices,
            const bool removeCells = false
        );


public:

    // Constructors

        //- Construct from polyMeshGen
        checkIrregularSurfaceConnections(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        checkIrregularSurfaceConnections
        (
            const checkIrregularSurfaceConnections&
        ) = delete;


    //- Destructor
    ~checkIrregularSurfaceConnections();

    // Member Functions

        //- Find vertices where the surface is not valid
        void checkIrregularVertices(labelHashSet& badVertices);

        //- Find invalid connections and fix them on the fly
        bool checkAndFixIrregularConnections();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const checkIrregularSurfaceConnections&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
