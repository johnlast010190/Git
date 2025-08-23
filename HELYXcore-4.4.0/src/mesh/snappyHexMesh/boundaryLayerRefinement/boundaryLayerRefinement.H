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
    (c) 2023 Engys Ltd.

Application
    Foam::boundaryLayerRefinement

Description
    Add 2D refinement of boundary layer cells

\*---------------------------------------------------------------------------*/

#ifndef boundaryLayerRefinement_H
#define boundaryLayerRefinement_H

#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class polyTopoChangeMap;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class boundaryLayerRefinement Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayerRefinement
{
public:

    // Public data

private:

    // Private data

        //- primitive patch to refine cells
        const indirectPrimitivePatch& pp_;

        //- Mesh+surface
        meshRefinement& meshRefiner_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("boundaryLayerRefinement");

    // Constructors

        //- Construct from mesh.
        boundaryLayerRefinement
        (
            const indirectPrimitivePatch& pp,
            meshRefinement& meshRefiner
        );

    // Member Functions

        // Helper

        void addTopFaces
        (
            const fvMesh& mesh,
            const labelList& pointLevel,
            const labelList& newEdgePts,
            const label& cLevel,
            const face& f,
            const label& anchorPt,
            const label& newMidPt,
            DynamicList<label>& faceVerts
        );

        void addAnchorToMid
        (
            const fvMesh& mesh,
            const labelList& pointLevel,
            const labelList& newEdgePts,
            const label& cLevel,
            const face& f,
            const label& anchorPt,
            DynamicList<label>& faceVerts
        );

        label findMidPoint
        (
            const fvMesh& mesh,
            const labelList& pointLevel,
            const labelList& newEdgePts,
            const label& cLevel,
            const face& f,
            const label& anchorStart,
            const label& anchorEnd
        );

        // Access

        // refine boundary cells

        void setRefinement();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
