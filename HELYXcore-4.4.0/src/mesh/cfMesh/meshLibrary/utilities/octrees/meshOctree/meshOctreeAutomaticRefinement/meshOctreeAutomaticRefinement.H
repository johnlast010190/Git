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
    meshOctreeAutomaticRefinement

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAutomaticRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAutomaticRefinement_H
#define meshOctreeAutomaticRefinement_H

#include "utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifier.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurfaceCurvatureEstimator;
class triSurfacePartitioner;
class meshOctreeAddressing;

/*---------------------------------------------------------------------------*\
                Class meshOctreeAutomaticRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAutomaticRefinement
{
    // Private data

        //- Reference to meshOctree
        meshOctree& octree_;

        //- Const reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- Use DATA boxes
        bool useDATABoxes_;

        //- A flag for activating paired refinement
        bool hexRefinement_;


    // Private member functions

        //- MeshOctreeAddressing pointer
        mutable meshOctreeAddressing* octreeAddressingPtr_;

        //- Calculate octree addresing
        void createOctreeAddressing() const;
        const meshOctreeAddressing& octreeAddressing() const;

        //- Surface partitioner pointer
        mutable triSurfacePartitioner* partitionerPtr_;

        //- Calculate partitionerPtr_
        void createSurfacePartitioner() const;
        const triSurfacePartitioner& partitioner() const;

        //- Curvature estimator pointer
        mutable triSurfaceCurvatureEstimator* curvaturePtr_;

        //- Calculate curvaturePtr_
        void createCurvatureEstimator() const;
        const triSurfaceCurvatureEstimator& curvature() const;


    // Functions and data performing automatic refinement

        //- Maximum allowed refinement level
        direction maxRefLevel_;

        //- Set maxRefLevel_
        void setMaxRefLevel();

        //- Refine selected boxes
        void refineSelectedBoxes(labelList&, labelLongList&);

        //- Refine DATA boxes based on curvature
        bool curvatureRefinement();
        bool refineBasedOnCurvature(labelList&, const labelLongList&);

        //- Refine DATA boxes which contain two or more independent parts
        //  of the surface in their vicinity
        bool proximityRefinement();
        bool refineBasedOnProximityTests
        (
            labelList&,
            const labelLongList&
        );

        //- Refine DATA boxes containing 2 or more corners, 2 or more distinct
        //  edges, etc.
        bool refineBasedOnContainedPartitions
        (
            labelList&,
            const labelLongList&
        );

        //- Refine boxes based on the number of contained surface corners
        bool refineBasedOnContainedCorners
        (
            labelList&,
            const labelLongList&
        );


public:

    // Constructors

        //- Construct from meshOctree, dictionary and boolean
        meshOctreeAutomaticRefinement
        (
            meshOctree& mo,
            const IOdictionary& dict,
            bool useDATABoxes = false
        );

        //- Disallow default bitwise copy construction
        meshOctreeAutomaticRefinement
        (
            const meshOctreeAutomaticRefinement&
        ) = delete;


    //- Destructor
    ~meshOctreeAutomaticRefinement();


    // Member Functions

        //- Activate hex refinement
        void activateHexRefinement();

        //- Refine octree boxes which are too coarse and cannot capture
        //  surface features
        void automaticRefinement();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAutomaticRefinement&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
