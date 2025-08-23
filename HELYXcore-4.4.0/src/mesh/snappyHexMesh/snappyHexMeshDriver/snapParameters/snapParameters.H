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
    (c) 2010-2024 Engys Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::snapParameters

Description
    Simple container to keep together snap specific information.

SourceFiles
    snapParameters.C

\*---------------------------------------------------------------------------*/

#ifndef snapParameters_H
#define snapParameters_H

#include "db/dictionary/dictionary.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/bools/lists/boolList.H"
#include "primitives/ints/lists/labelList.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyBoundaryMesh;
class refinementSurfaces;
/*---------------------------------------------------------------------------*\
                       Class snapParameters Declaration
\*---------------------------------------------------------------------------*/

class snapParameters
{
    // Private data

        const label nSmoothPatch_;

        const label nSmoothInternal_;

        scalar snapTol_;

        const label nSmoothDispl_;

        const label nSnap_;

        //- Per patch whether to perform feature snapping.
        boolList featureEdgePatches_;

        //- Per patch whether to perform region snapping.
        boolList regionSnapPatches_;

        //- Per patch whether to perform direct snapping.
        boolList directSnapPatches_;

        // Enable/disable direct feature snapping
        bool directFeatureSnapping_;

        // Repartition before direct feature snapping
        bool directSnapRepartition_;

        // Enable/disable feature snapping for zoned surfaces
        bool zoneFeatureSnapping_;

        // Setting for correcting concave faces (-1 to disable)
        const scalar concaveTol_;

        // Setting for edge collapsing (-1 to disable)
        const scalar collapseTol_;

        // Enable pre feature snapping iterations
        label nPreFeatureIter_;

        // Enable main feature snapping iterations
        label nFeatureIter_;

        // Number of outer feature snapping loops
        label nOuterIter_;

        bool globalFE_;

        bool globalDirectSnap_;

        bool globalSnapRegions_;

        bool globalPreMerge_;

        //Split boundary cells where face shared by more than one edge
        bool splitDegenerateCells_;

        const Switch explicitFeatureSnap_;

        const Switch implicitFeatureSnap_;

        const Switch multiRegionFeatureSnap_;

        const Switch detectNearSurfacesSnap_;

        const Switch strictRegionSnap_;

        const Switch detectBaffles_;

        const Switch baffleFeaturePoints_;

        const Switch releasePoints_;

        const Switch stringFeatures_;

        const Switch avoidDiagonal_;

        //- How often needs face splitting be run
        label nFaceSplitInterval_;

        //- When is angle too concave too split
        scalar concaveAngle_;

        //- When is face-split not sufficiently diagonal
        scalar minAreaRatio_;

        //Compute more expensive surface normal calculation during projection
        const Switch averageSurfaceNormal_;

        //Switch to enable larger stencil for surface averaging
        const Switch enlargeStencil_;

        //Switch to enable smoothing of surface during feature snapping
        const Switch smoothSnappedSurface_;

        //whether to re-order baffles for generating cyclics
        const Switch reorderBaffles_;

        //Point-face smoothing of displacement during feature snapping
        const label nSmoothFeatureDisp_;

        //Number of smoothing iteration for sliver faces
        const label nSliverSmooths_;

        //Try adding additional tetrahedral cells to split mesh
        const Switch addTetsToSplitMesh_;

        //Perform error checks during feature snap
        const Switch featureSnapChecks_;

        //Write vtk files of snapping process
        const Switch writeSnapVTK_;

        //Number of additional layers to add after creating dual mesh
        const label nAdditionalDualLayers_;

        //Perform initial snap to zoned surfaces
        const Switch preZoneSnap_;

        //Perform merging of boundary faces where possible
        const Switch mergeBoundaryFaces_;

        //Perform merging of boundary faces if they are on different patches
        const Switch mergeAcrossPatches_;

        //Repatch overlapping surfaces to have lowest ID
        const Switch repatchOverlapping_;

        //Added pre-smoothing of nodes adjacent to baffle edges
        const Switch preSmoothBaffles_;

        //Additional feature snapping for acute and reflex angles
        const scalar acuteReflexSnap_;

        //Additional feature snapping for acute and reflex angles (minimum)
        const scalar minAcuteReflexSnap_;

        //Measure of conformity to surface relative to mesh spacing
        // used for smoothing non-conformal points
        const scalar conformityTol_;

        //For poor conformity the number of smoothing iterations
        const label nConformitySmooths_;

        //Check how well boundary mesh conforms to the surface
        const Switch checkConformity_;

        //Check whether to use old dual layer removal method
        const Switch dualLayerRemoval_;

        //Whether to lock onto weak feature angles in specified range (degrees)
        const Tuple2<scalar,scalar> weakFeatureSnap_;

        //Optional merging boundary faces before extruding
        const Switch preMergeExtrude_;

        //Optional per surface testing of nearest snap performance
        const Switch testPerformance_;

        //Optional iteration to perform flipping of non-convergent faces
        const label snapFlipIter_;

        //Optional tolerance used in deciding to flip faces
        const scalar snapFlipTol_;

        //Optional merging faces on grown up patches
        const Switch mergeGrownUp_;

public:

    // Constructors

        //- Construct from dictionary - new syntax with additional features
        snapParameters
        (
            const dictionary& dict,
            const polyBoundaryMesh&
        );

        //- Disallow default bitwise copy construction
        snapParameters(const snapParameters&) = delete;


    // Member Functions

            void calculateRegionAndFeatureSnapPatches
            (
                const dictionary& dict,
                const meshRefinement& meshRefiner
            );

        // Access

            //- Number of patch smoothing iterations before finding
            //  correspondence to surface
            label nSmoothPatch() const
            {
                return nSmoothPatch_;
            }

            //- Number of internal point smoothing iterations (combined with
            //  nSmoothPatch
            label nSmoothInternal() const
            {
                return nSmoothInternal_;
            }

            //- Relative distance for points to be attracted by surface
            //  feature point
            //  or edge. True distance is this factor times local
            //  maximum edge length.
            scalar snapTol() const
            {
                return snapTol_;
            }
            scalar& snapTol()
            {
                return snapTol_;
            }

            //- Number of mesh displacement smoothing iterations.
            label nSmoothDispl() const
            {
                return nSmoothDispl_;
            }

            //- Maximum number of snapping relaxation iterations. Should stop
            //  before upon reaching a correct mesh.
            label nSnap() const
            {
                return nSnap_;
            }

            // Surface snapping specific

                //- Override attraction to nearest with intersection location
                //  at near surfaces
                Switch detectNearSurfacesSnap() const
                {
                    return detectNearSurfacesSnap_;
                }

                //- Attract point to corresponding surface region only
                Switch strictRegionSnap() const
                {
                    return strictRegionSnap_;
                }


            // Feature edge snapping specific

                //- Maximum number of pre feature snapping relaxation iterations.
                label nPreFeatureIter() const
                {
                    return nPreFeatureIter_;
                }

                Switch explicitFeatureSnap() const
                {
                    return explicitFeatureSnap_;
                }

                Switch implicitFeatureSnap() const
                {
                    return implicitFeatureSnap_;
                }

                Switch multiRegionFeatureSnap() const
                {
                    return multiRegionFeatureSnap_;
                }

                Switch detectBaffles() const
                {
                    return detectBaffles_;
                }

                Switch baffleFeaturePoints() const
                {
                    return baffleFeaturePoints_;
                }

                Switch releasePoints() const
                {
                    return releasePoints_;
                }

                Switch stringFeatures() const
                {
                    return stringFeatures_;
                }

                Switch avoidDiagonal() const
                {
                    return avoidDiagonal_;
                }


                // Face splitting

                    label nFaceSplitInterval() const
                    {
                        return nFaceSplitInterval_;
                    }

                    scalar concaveAngle() const
                    {
                        return concaveAngle_;
                    }

                    scalar minAreaRatio() const
                    {
                        return minAreaRatio_;
                    }




            Switch enlargeStencil() const
            {
                return enlargeStencil_;
            }

            Switch smoothSnappedSurface() const
            {
                return smoothSnappedSurface_;
            }

            Switch averageSurfaceNormal() const
            {
                return averageSurfaceNormal_;
            }

            Switch reorderBaffles() const
            {
                return reorderBaffles_;
            }

            label nSmoothFeatureDisp() const
            {
                return nSmoothFeatureDisp_;
            }

            label nSliverSmooths() const
            {
                return nSliverSmooths_;
            }

            Switch addTetsToSplitMesh() const
            {
                return addTetsToSplitMesh_;
            }

            Switch featureSnapChecks() const
            {
                return featureSnapChecks_;
            }

            Switch writeSnapVTK() const
            {
                return writeSnapVTK_;
            }

            label nAdditionalDualLayers() const
            {
                return nAdditionalDualLayers_;
            }

            Switch preZoneSnap() const
            {
                return preZoneSnap_;
            }

            Switch mergeBoundaryFaces() const
            {
                return mergeBoundaryFaces_;
            }

            Switch mergeAcrossPatches() const
            {
                return mergeAcrossPatches_;
            }

            Switch repatchOverlapping() const
            {
                return repatchOverlapping_;
            }

            Switch preSmoothBaffles() const
            {
                return preSmoothBaffles_;
            }

            scalar acuteReflexSnap() const
            {
                return acuteReflexSnap_;
            }

            scalar minAcuteReflexSnap() const
            {
                return minAcuteReflexSnap_;
            }

            scalar conformityTol() const
            {
                return conformityTol_;
            }

            label nConformitySmooths() const
            {
                return nConformitySmooths_;
            }

            Switch checkConformity() const
            {
                return checkConformity_;
            }

            Switch dualLayerRemoval() const
            {
                return dualLayerRemoval_;
            }

            Tuple2<scalar, scalar> weakFeatureSnap() const
            {
                return weakFeatureSnap_;
            }

            Switch preMergeExtrude() const
            {
                return preMergeExtrude_;
            }

            Switch testPerformance() const
            {
                return testPerformance_;
            }

            label snapFlipIter() const
            {
                return snapFlipIter_;
            }

            scalar snapFlipTol() const
            {
                return snapFlipTol_;
            }

            Switch mergeGrownUp() const
            {
                return mergeGrownUp_;
            }

            //- Maximum number of feature snapping relaxation iterations.
            label nFeatureIter() const
            {
                return nFeatureIter_;
            }

            //- Number of outer feature snapping loops.
            label nOuterIter() const
            {
                return nOuterIter_;
            }

            // Whether we want feature edges everywhere
            bool globalFE() const
            {
                return globalFE_;
            }

            // Whether we want to perform direct feature snapping everywhere
            bool globalDirectSnap() const
            {
                return globalDirectSnap_;
            }

            // Whether we want region snapping everywhere
            bool globalSnapRegions() const
            {
                return globalSnapRegions_;
            }

            // Whether we want to pre-snap everwhere (extrude algorithm)
            bool globalPreMerge() const
            {
                return globalPreMerge_;
            }

            // Whether to face split degenerate cells
            bool splitDegenerateCells() const
            {
                return splitDegenerateCells_;
            }

            // Whether we want to perform direct feature snapping
            bool directFeatureSnapping() const
            {
                return directFeatureSnapping_;
            }

            // Whether to repartition before perform direct feature snapping
            bool directSnapRepartition() const
            {
                return directSnapRepartition_;
            }

            // Whether we want to perform zonal feature snapping
            bool zoneFeatureSnapping() const
            {
                return zoneFeatureSnapping_;
            }

            // Tolerance used for concavity test
            scalar concaveTol() const
            {
                return concaveTol_;
            }

            // Tolerance used for edge collapsing
            scalar collapseTol() const
            {
                return collapseTol_;
            }

            labelList regionSnappedIDs() const
            {
                DynamicList<label> snapRegions(regionSnapPatches_.size());

                forAll(regionSnapPatches_, i)
                {
                    if (regionSnapPatches_[i])
                    {
                        snapRegions.append(i);
                    }
                }
                snapRegions.shrink();

                return labelList(snapRegions, true);
            }

            labelList directSnappedIDs() const
            {
                DynamicList<label> directRegions(directSnapPatches_.size());

                forAll(directSnapPatches_, i)
                {
                    if (regionSnapPatches_[i] || directSnapPatches_[i])
                    {
                        directRegions.append(i);
                    }
                }
                directRegions.shrink();

                return labelList(directRegions, true);
            }

        // Per patch information

            //- Per patch whether to perform feature edge snapping.
            const boolList& featureEdgePatches() const
            {
                return featureEdgePatches_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const snapParameters&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
