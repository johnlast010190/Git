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
    (c) 2015-2024 Engys Ltd

Class
    Foam::dynamicGIBFvMesh

Description
    Class for snapping internal point to an STL and making the faceZone needed
    for the GIB

SourceFiles
    dynamicGIBFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicGIBFvMesh_H
#define dynamicGIBFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "fields/volFields/volFields.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "primitives/Tuple2/Tuple2.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapGIB;

/*---------------------------------------------------------------------------*\
                           Class dynamicGIBFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicGIBFvMesh
:
    public dynamicFvMesh
{
    // Private data

    class xxorEqOp
    {
        public:

        void operator()(bool& x, const bool& y) const
        {
            if (!x && y)
            {
                x = false;
            }
            else if (x && !y)
            {
                x = false;
            }
            else
            {
                x = true;
            }
        }
    };

    //Operator copied from meshRefinement
    class maxMagEqOp
    {

    public:

        void operator()(vector& x, const vector& y) const
        {
            for (direction i = 0; i < vector::nComponents; i++)
            {
                scalar magX = mag(x[i]);
                scalar magY = mag(y[i]);

                if (magX < magY)
                {
                    x[i] = y[i];
                }
                else if (magX == magY)
                {
                    if (y[i] > x[i])
                    {
                        x[i] = y[i];
                    }
                }
            }
        }
    };


    //- Private functions

        //- Initialise after construction
        virtual void initialize();

        virtual void modifyRegion0Patch();

        //- Create GIB face addressing
        virtual void makeFl() const;

        //- Create GIB flipmap
        virtual void makeFlipMap() const;

        virtual void makeInterPoints() const;

        virtual void makeInterPoints0() const;

        virtual void makePhiB() const;

        virtual void makePhiTr() const;

        virtual void makePhiTrS() const;

        virtual void makePhiTrG() const;

        virtual void makeFaceIndicator() const;

        virtual void makeFaceIndicator0() const;

        virtual void makecRegion() const;

        virtual void makecRegion0() const;

        virtual void makeFullSnapCells() const;

        virtual void makeFullSnapCellPoints() const;

        virtual void makePopCellPoints() const;

        virtual void makePopSPoints() const;

        virtual void makePopGPoints() const;

        virtual void makeConcavePoints() const;

        virtual void makeCloseBoundaryPoints() const;

        virtual void makeMarkedBoundaryPoints() const;

        virtual void makeNormal2DEdges() const;

        virtual void makeBoundaryPointsRev() const;

        //- reconstruct old point location via stl mapping
        virtual void makeRecPoints0() const;

        virtual void makeRecAllPoints0() const;

        virtual void makePopUpCells() const;

        virtual void makeBoundaryPoints() const;

        void extendBoundaryBaseMesh(const labelList&, vectorField&) const;


protected:


    // Protected Member variables

        //- Dictionary of motion control parameters
        const dictionary dynamicMeshCoeffs_;

        //- The reference points which are transformed
        mutable pointIOField* basePoints_;

        //- The reference points which are transformed
        mutable vectorField* baseCC_;

        //- The reference points which are transformed
        mutable vectorField* baseSf_;

        //- The reference points which are transformed
        mutable vectorField* baseCf_;

        //- Name of the stl geometry
        const word triName_;

        //- GIB patches
        mutable boolList* patchGIBPtr_;

        //- Master GIB patches
        label masterGIB_;

        //- Slave GIB patches
        label slaveGIB_;

        //- Mesh regions label
        mutable labelList* cRegionPtr_;

        //- Mesh regions label
        mutable labelList* cRegion0Ptr_;

        //- Points connected to boundary faces via edges
        mutable labelListList* closeBoundaryPointsPtr_;

        //- Points connected to boundary faces via edges
        mutable boolList* markedBoundaryPointsPtr_;

        //- Edges normal to the 2d planes
        mutable boolList* normal2DEdgesPtr_;

        //- Points to be checked for smnapping to base location
        mutable boolList* boundaryPointsRevPtr_;

        //- Stl geometry
        mutable triSurfaceMesh* ibMeshPtr_;

        //- Interface faceZone addressing
        mutable labelList* flPtr_;

        //- Interface faceZone flipmap
        mutable boolList* fmPtr_;

        //- True/false for snapping points
        mutable boolList* interPointsPtr_;

        //- True/false for snapping points
        mutable boolList* interPoints0Ptr_;

        //- True/false for snapping points
        mutable boolList* boundaryPointsPtr_;

        //- Faces with more than 1 intersection
        mutable boolList* multInterFacesPtr_;

        //- Mark the concave points
        mutable boolList* concavePointsPtr_;

        //- Mark the boundary edges around a concave point
        mutable boolList* concaveEdgesPtr_;

        //- List of cells that have all the points are snapp points
        mutable labelList* fullSnapCellsPtr_;

        //- Points of cells that have all the points snapped
        mutable labelList* fullSnapCellPointsPtr_;

        //- Identifier (global) for points of pop cells
        mutable boolList* popCellPointsPtr_;

        //- Identifier (global) for points that snap during shrink stage
        mutable boolList* popSPointsPtr_;

        //- Identifier (global) for points that unsnap during grow stage
        mutable boolList* popGPointsPtr_;

        //- Triangles of stl
        //  needed only from transformation class
        mutable labelIOList* hitIndexPtr_;
        mutable vectorField* hitPointPtr_;

        //- Identifier(global) for pop cells
        mutable boolList* popUpCellsPtr_;

        //- Identifier (global) for faces which are part of the interface
        mutable boolList* faceIndicatorPtr_;

        //- Identifier (global) for faces which are part of the interface
        //  in the previous timestep
        mutable boolList* faceIndicator0Ptr_;

        //- Visualization - faceCells of side one
        boolList masterFCells_;

        //- Visualization - faceCells of side two
        boolList slaveFCells_;

        //- Patch names that "touch" the fluid/region0
        wordList region0Patch_;

        //- Fields that limit if negative
        wordList limitNegFieldNames_;

        //- Number of popping sub-steps
        label popSubSteps_;

        //- The maximum allowed non-orthogonality stored as cosine of the angle
        scalar maxNonOrth_;

        //- Unsnap percentage of fully snapped cells
        scalar unsnapVar_;

        //- Old time faceZone addressing
        mutable labelList* fl0Ptr_;

        //- Old time faceZone flipmap
        mutable boolList* fm0Ptr_;

        //- Mesh flux due to boundary motion
        mutable surfaceScalarField* phiBPtr_;

        //- Total mesh flux due to transition
        mutable surfaceScalarField* phiTrPtr_;

        //- Mesh flux due to shring stage
        mutable surfaceScalarField* phiTrSPtr_;

        //- Mesh flux due to grow stage
        mutable surfaceScalarField* phiTrGPtr_;

        //- Switch for debug mode
        bool debugMode_;

        //-
        bool correctBCs_;

        bool pyrPrismFlip_;

        //- Clip pop cell values with min/max of surrounding faces
        bool boundPopValues_;

        //- Check if interface has the same components
        // as previous timesteps
        bool checkInterface_;

        //- Reconstructed old points
        pointField oldPoints_;

        //- Point location of previous timestep
        pointField prevPoints_;

        //- Reconstructed old points
        mutable pointField* recPoints0Ptr_;

        mutable pointField* recAllPoints0Ptr_;

        //- Indicator for points which are moved due to quality controls
        mutable boolList* meshQualityPointsPtr_;


    // Protected Member Functions

        virtual bool includeWalls() const = 0;

        //- Do tangential Boundary Motion
        virtual bool tangentialBoundaryMotion() const = 0;


        //- Invert patch to have one side normals
        virtual faceList preparePatch(const fvPatch&);

        void syncPoints(pointField&) const;

        virtual labelList checkingIntersectingFaces() const;

        virtual bool checkIfEdgeToleranceIntersecting
        (
            const face&,
            const face&
        ) const;

        virtual void doTangentialBoundaryMotion
        (
            pointField&
        ) const;


        virtual tmp<pointField> findSnappedPoints
        (
            const bool& fromBase = true
        ) const;

        //- Roll back the snapping as necessary to satisfy quality constraint
        virtual void unsnapForQuality
        (
            pointField&,
            const pointField&
        ) const;

        virtual void fullSnappedPointsTreatment
        (
            pointField&,
            const pointField&
        ) const;

        //- Approximate old positions of unsnapped cells points
        virtual void computeOldPositionsInUnsnappedCells
        (
            pointField&,
            const labelList&
        ) const;


        //- Returns the label of the faceZone.
        //  if doesn't find (index==-1), it creates it
        virtual label gibFaceZone() const;

        //- Returns the label of the cellZone.
        //  if doesn't find (index==-1), it creates it
        virtual label findOrCreateCellZone() const;

        //- If boundary points (which are on GIB) are moved tangentialy
        //  to the boundary face and they dont have an internal edge
        //  moving to their original location
        virtual void correctBoundaryPointsOnBaseMesh
        (
            pointField&
        );



        //- Algorithm for creating the flipmap faster
        virtual void calculateFlipMap() const;

        //- Algorithm for creating the flipmap faster
        virtual void regionVisDebug() const;

        virtual boolList findMultInterPoints() const;

        //- Distribute fields - Shrink stage
        virtual void popShrinkFields();

        //- Distribute fields - Grow stage
        virtual void popGrowFields();

        //- Compute the old volumes based on the current volume
        //  and the mesh fluxes
        virtual void correctV0();

        //- Function takes small/negative V0 volumes that were
        //  corrected to small and corrects fluxes to correspond
        //  to V0 correctly
        virtual void negativeV0CorrectGIBFluxes
        (
            surfaceScalarField& phiB,
            const scalarField& V,
            const scalarField& V0,
            const dimensionedScalar& deltaT
        );

        //- Modify old, current volumes and mesh fluxes of squeezed cells
        //  Squeezed are cells with vol < SMALL
        virtual void modifyVolSqueezedCells() const;

        virtual void resetMeshFluxes();
        virtual void correctVelocityFlux();

        //- Move old GIB
        virtual tmp<vectorField> movePolyPatch
        (
            primitivePatch&,
            const vectorField&
        ) = 0;


        virtual void faceCellsVisDebug();
        virtual void visCells();
        virtual void findGIBPatches();
        virtual void findGIBPatches(const word&);
        virtual void modifyRegionLabels(labelList&) const;

        virtual bool isIdenticalInterface() const;

        //- Smooth castelated surface to improve snapping
        virtual void smoothInternalBasePoints
        (
            pointField&,
            const labelList&
        ) const;


        // Functions for flipping cells

            //- Flipcells true/false (solid/fluid side)
            void flipCells(boolList&, bool) const;

            void singleCellSuddenPopCellRemoval(boolList&) const;

            //- Find all cells with all their points at the interface
            boolList checkFlipCell(const boolList&) const;

            //- 0 --> nothing
            //- 1 --> now not interface | change to interface
            //- -1 --> now interface | change to non-interface
            //- 2 between flipping cells | now not interface
            //- -2 between flipping cells | now interface
            labelList faceFlippingType
            (
                const boolList&,
                const boolList&
            ) const;


            labelList cellFlippingNumber
            (
                const boolList&,
                const boolList&
            ) const;


            void findFlippingFaces
            (
                const boolList&,
                boolList&,
                const boolList&,
                const labelList&,
                labelList&
            ) const;

            void resetFlipCellsWithFewInterfaces
            (
                boolList&,
                const labelList&
            ) const;

            void resetFlipHexCellsWith3Interfaces
            (
                boolList&,
                const labelList&,
                labelList&
            ) const;


        // Functions only for debugging

            //- Writes the cells with face info
            virtual void writeProblematicCells(label) const;

            //- Writes the cells with face info + field
            virtual void writeProblematicCells
            (
                label,
                const surfaceScalarField&
            ) const;

            //- Writes all the cells that have either zero V or V0
            virtual void writeProblematicCells() const;

            //- Quality control
            virtual void qualityControl() const;

            //- Writes all the cells that have all points snapped
            virtual void writeProblematicCells(const boolList&) const;

            //- Check if didn't find different regions
            virtual void checkRegion(const labelList&) const;

            //- Check if didn't find different regions
            virtual void writeRegionInterface
            (
                word,
                const boolList&
            ) const;

            //- Check if didn't find different regions
            virtual void writeRegionInterface
            (
                word,
                const boolList&,
                label
            ) const;

            //- Provide a scalarfield and write it as an volScalarField
            //  for post processing
            template<class Type>
            void writeScalarField
            (
                word,
                const List<Type>&,
                bool writeNow=false
            ) const;

        virtual void checkConcaveBoundaryPoints
        (
            const pointField&,
            pointField&
        ) const;

        virtual void correctConstraintPatches
        (
            pointField&
        ) const;

        virtual void correctSymmetryPatches
        (
            pointField&
        ) const;

        virtual void applyTwoDCorrection
        (
            pointField&
        )const;

        virtual void applyTwoDPlanesCorrection
        (
            pointField&,
            pointField&
        );

        //- Create inactive cell Zone
        virtual void updateSolidCellZone() const;

        //- Store the old time fields
        virtual void storeOldTimes();

        virtual void correctBCs();

        virtual void explicitLimitingFields();

        //- Clear all demand-driven data
        virtual void clearOutGIBData();

        //- Shrink fields pop-cells treatment
        template<class Type, template<class> class PatchField, class GeoMesh>
        void shrinkPopFields
        (
            const labelList& popShrinkFaces,
            const surfaceScalarField& popShrinkPhi,
            const boolList& popShrinkCells
        );

        //- Grow fields pop-cells treatment
        template<class Type, template<class> class PatchField, class GeoMesh>
        void growPopFields
        (
            const labelList& popGrowFaces,
            const surfaceScalarField& popGrowPhi,
            const boolList& popGrowCells,
            const scalarField& oldVolgrow
        );


public:

    //- Runtime type information
    TypeName("dynamicGIBFvMesh");


        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicGIBFvMesh,
            IOobject,
            (
                const IOobject& io,
                const dictionary dict
            ),
            (io, dict)
        );


    // Constructors

        //- Construct from IOobject
        dynamicGIBFvMesh(const IOobject&, const word&);

        //- Construct from IOobject
        dynamicGIBFvMesh
        (
            const IOobject&,
            const dictionary
        );

        //- Disallow default bitwise copy construction
        dynamicGIBFvMesh(const dynamicGIBFvMesh&) = delete;


    //- Destructor
    ~dynamicGIBFvMesh();


    // Member Functions

        //-
        virtual const triSurfaceMesh& triSM() const
        {
            if (!ibMeshPtr_)
            {
                FatalErrorInFunction
                    << "Internal boundary surface mesh not set" << nl
                    << exit(FatalError);
            }
            return *ibMeshPtr_;
        }

        label masterId() const
        {
            return masterGIB_;
        }

        label slaveId() const
        {
            return slaveGIB_;
        }


        //- FaceZone initialization using the stl information
        virtual void updateInit(const word&);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() override
        {
            return false;
        }

        //- Encapsulated update algorithm
        virtual void doUpdate(mapGIB&, bool);

        virtual void writeGeometry(const fileName& surfaceName);

        //- Correct boundaryConditions
        virtual tmp<vectorField> velocityCorrect(const vectorField&) const override = 0;

        //- Calculate the old boundary location
        virtual tmp<vectorField> oldBoundaryLocation() const = 0;

        //-
        virtual const boolList& interP() const;

        //-
        virtual const boolList& interP0() const;

        //-
        virtual const labelList& fl() const;

        //-
        virtual const labelList& fl0() const;

        //-
        virtual const boolList& fm() const;

        virtual const boolList& fm0() const;

        //- Mesh flux due to boundary motion
        virtual const surfaceScalarField& phiB() const;

        //- Total mesh flux due to transition
        virtual const surfaceScalarField& phiTr() const;

        //- Mesh flux due to shring stage
        virtual const surfaceScalarField& phiTrS() const;

        //- Mesh flux due to grow stage
        virtual const surfaceScalarField& phiTrG() const;

        virtual const boolList& faceIndicator() const;

        virtual const boolList& faceIndicator0() const;

        //- Mesh regions label (ususually labels inactivated cells)
        virtual const labelList& cRegion() const;

        //- Mesh regions label at time 0 (ususually labels inactivated cells)
        virtual const labelList& cRegion0() const;

        virtual const labelList& fullSnapCells() const;

        virtual const labelList& fullSnapCellPoints() const;

        virtual const boolList& popCellPoints() const;

        virtual const boolList& popSPoints() const;

        virtual const boolList& popGPoints() const;

        virtual const boolList& concavePoints() const;

        virtual const boolList& concaveEdges() const;

        virtual const labelListList& closeBoundaryPoints() const;

        virtual const boolList& markedBoundaryPoints() const;

        virtual const boolList& boundaryPointsRev() const;

        virtual const pointField& recPoints0() const;

        virtual const pointField& recAllPoints0() const;

        virtual const boolList& popUpCells() const;

        virtual const boolList& boundaryPoints() const;

        virtual const labelIOList& hitIndex() const
        {
            return *hitIndexPtr_;
        }

        virtual const vectorField& hitPoint() const
        {
            return *hitPointPtr_;
        }

        //- Acess to edges normal to the 2d planes
        virtual const boolList& normal2DEdges() const
        {
            return *normal2DEdgesPtr_;
        }

        //- Acess to indicator for points which are moved
        //  due to quality controls
        virtual const boolList& meshQualityPoints() const
        {
            return *meshQualityPointsPtr_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicGIBFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
