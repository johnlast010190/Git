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
    (c) 2021 Engys Ltd.

Application
    Foam::autoFanExtrude

Description
    Extrude fan mesh from set of internal faces

\*---------------------------------------------------------------------------*/

#ifndef autoFanExtrude_H
#define autoFanExtrude_H

#include "meshRefinement/meshRefinement.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class polyTopoChangeMap;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class autoFanExtrude Declaration
\*---------------------------------------------------------------------------*/

class autoFanExtrude
{
        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        meshDistribute& distributor_;

        //- Container for layer information
        const layerParameters& layerParams_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("autoFanExtrude");

    // Constructors

        //- Construct from mesh.
        autoFanExtrude
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            meshDistribute& distributor,
            const layerParameters& layerParams
        );

    // Member Functions

        // Helper

        // Access

        // extrude from set of internal faces

            //- Perform fan extrusion
            void setRefinement();

            //- Identify faces to extrude
            bool facesToExtrude
            (
                const scalar& fanAngleCos,
                boolList& visited
            );

            //- Rebalance so faces to baffle on single processor
            void distribute
            (
                boolList& visited
            );

            //- Filter out unsuitable faces and split corner faces
            bool filterAndSplit
            (
                boolList& visited
            );

            //- create baffle faces
            void createBaffles
            (
                const boolList& visited,
                List<labelPair>& baffles
            );

            //- Duplicate baffle points
            void duplicate
            (
                labelList& bafflePtType,
                scalarField& minExternalLen,
                List<labelPair>& baffles
            );

            //- move points to create void for layer mesh
            void move
            (
                const indirectPrimitivePatch& bafflePP,
                const labelList& bafflePtType,
                const scalarField& minExternalLen,
                vectorField& disp
            );

            //- Extrude single layer mesh
            autoPtr<polyTopoChangeMap> extrudeLayer
            (
                const indirectPrimitivePatch& bafflePP,
                const labelList& bafflePtType,
                const vectorField& disp
            );

            //- Merge baffle faces
            void mergeBaffles
            (
                const polyTopoChangeMap& extrudeMap,
                List<labelPair>& baffles
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
