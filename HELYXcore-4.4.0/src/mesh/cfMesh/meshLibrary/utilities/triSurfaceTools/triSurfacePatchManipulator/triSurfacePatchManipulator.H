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
    triSurfacePatchManipulator

Description
    Generates patches in the surface mesh based
    on the user-selected feature edges

SourceFiles
    triSurfacePatchManipulator.C
    triSurfacePatchManipulatorFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfacePatchManipulator_H
#define triSurfacePatchManipulator_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfacePatchManipulator Declaration
\*---------------------------------------------------------------------------*/

class triSurfacePatchManipulator
{
    // Private data

        //- Const reference to triSurf
        const triSurf& surf_;

        //- Detected feature edges
        List<direction> featureEdges_;

        //- Surface patches
        labelList facetInPatch_;

        //- Number of patches
        label nPatches_;

        //- Patch names
        wordList newPatchNames_;

        //- Patch types
        wordList newPatchTypes_;


    // Private member functions

        //- Allocate and fill the feature edges list
        void allocateFeatureEdges();

        //- Create patches bounded by a set of feature edges
        void createPatches();


public:

    // Constructors

        //- Construct from triSurface
        triSurfacePatchManipulator(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfacePatchManipulator(const triSurfacePatchManipulator&) = delete;


    //- Destructor
    ~triSurfacePatchManipulator();


    // Member Functions

        void detectedSurfaceRegions(VRWGraph&) const;

        //- Store regions into subsets with a given prefix
        const triSurf* surfaceWithPatches
        (
            IOdictionary* meshDictPtr = nullptr,
            const word prefix = "patch_",
            const bool forceOverwrite = false
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfacePatchManipulator&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
