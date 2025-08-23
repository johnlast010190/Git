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
    (c) 2015-2024 Engys Ltd.

Class
    Foam::solidBodyGIBFvMesh

Description
    Class for snapping internal point to an STL and making the faceZone needed
    for the GIB

SourceFiles
    solidBodyGIBFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyGIBFvMesh_H
#define solidBodyGIBFvMesh_H

#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "referenceFrames/coordinateFrame.H"
#include "GIB/geometryTransformation/geometryTransformation.H"
#include "fields/volFields/volFields.H"
#include "primitives/Tuple2/Tuple2.H"
#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidBodyGIBFvMesh Declaration
\*---------------------------------------------------------------------------*/

class solidBodyGIBFvMesh
:
    public dynamicGIBFvMesh
{
    // Private Member Data

        //- The previous time stl data (stl goes through multiple transforms)
        //  these points are just to allow incremental motion definition.
        mutable autoPtr<pointField> geomOldPoints_;


    // Private Member Functions

        void initialise();

        //- Compute the stl <--> GIB addressing
        //  Needed in case when GIB initialized using existing faceZone which is
        //  being triangulated via the utility
        void calculateHitIndex();


protected:

    // Protected data

        //- init stl points
        mutable pointField* initStlPointsPtr_;

        //- Is the motion solid body motion?
        bool isSolidBody_;

        //- The solid body motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Pointer to moving reference frame
        coordinateFrame* coorFramePtr_;

        //- The deformation of the surface
        PtrList<geometryTransformation> geoTrans_;

        autoPtr<motionSolver> motionPtr_;


    // Protected Member Functions

        virtual bool includeWalls() const
        {
            return false;
        }

        //- do tangential Boundary Motion
        virtual bool tangentialBoundaryMotion() const
        {
            return false;
        }

        //- move old GIB
        virtual tmp<vectorField> movePolyPatch
        (
            primitivePatch&,
            const vectorField&
        );

        virtual tmp<pointField> moveSurface(bool updateSolver=true);

        //- reconstruction of the new patch based on the surface movement
        virtual tmp<pointField> recNewPointLocation
        (
            const pointField&,
            const labelList&
        );

        virtual void computeOldPositionsInUnsnappedCells
        (
            pointField&,
            const labelList&
        ) const;

public:

    //- Runtime type information
    TypeName("solidBodyGIBFvMesh");


    // Constructors

        //- Construct from IOobject
        solidBodyGIBFvMesh(const IOobject&, const word&);

        //- Construct from IOobject
        solidBodyGIBFvMesh(const IOobject& io);

        //- Construct from IOobject
        solidBodyGIBFvMesh(const IOobject& io, const dictionary);

        //- Disallow default bitwise copy construction
        solidBodyGIBFvMesh(const solidBodyGIBFvMesh&) = delete;


    //- Destructor
    ~solidBodyGIBFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- faceZone initialization using the stl information
        virtual void updateInit(const word&);

        //- Correct boundaryConditions
        virtual tmp<vectorField> velocityCorrect(const vectorField&) const;

        //- Calculate the old boundary location
        virtual tmp<vectorField> oldBoundaryLocation() const;

        virtual const pointField& oldPoints() const
        {
            return oldPoints_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidBodyGIBFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
