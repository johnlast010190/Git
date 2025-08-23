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
    (c) 2024 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::solidBody

Description
    Dynamic GIB mesh changer for solid body motion.

SourceFiles
    fvMeshGIBChangersSolidBody.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshGIBChangersSolidBody_H
#define fvMeshGIBChangersSolidBody_H

#include "base/fvMeshGIBChangersBase.H"
#include "GIB/geometryTransformation/geometryTransformation.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                 Class fvMeshGIBChangers::solidBody Declaration
\*---------------------------------------------------------------------------*/

class solidBody
:
    public fvMeshGIBChangersBase
{
    // Private Member Data

        //- The previous time stl data (stl goes through multiple transforms)
        //  these points are just to allow incremental motion definition.
        mutable autoPtr<pointField> geomOldPoints_;

    // Private Member Functions

        void initialise();

        //- Compute the stl <--> GIB addressing.
        //  Needed in case when GIB initialized using existing faceZone,
        //  which is being triangulated via the utility.
        void calculateHitIndex();


protected:

    // Protected Data

        //- Init stl points
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

        virtual bool includeWalls() const override
        {
            return false;
        }

        //- Do tangential boundary motion
        virtual bool tangentialBoundaryMotion() const override
        {
            return false;
        }

        //- Move old GIB
        virtual tmp<vectorField> movePolyPatch
        (
            primitivePatch&,
            const vectorField&
        ) override;

        virtual tmp<pointField> moveSurface(bool updateSolver = true);

        //- Reconstruction of the new patch based on the surface movement
        virtual tmp<pointField> recNewPointLocation
        (
            const pointField&,
            const labelList&
        );

        virtual void computeOldPositionsInUnsnappedCells
        (
            pointField&,
            const labelList&
        ) const override;


public:

    //- Runtime type information
    TypeName("solidBody");


    // Constructors

        //- Construct from fvMesh
        explicit solidBody(fvMesh&);

        //- Construct from fvMesh and type
        solidBody(fvMesh&, const word&);

        //- Construct from fvMesh and dictionary
        solidBody(fvMesh&, const dictionary);

        //- Disallow default bitwise copy construction
        solidBody(const solidBody&) = delete;


    //- Destructor
    ~solidBody();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() override;

        //- faceZone initialization using the stl information
        virtual void updateInit(const word&) override;

        //- Correct velocity for GIB
        virtual tmp<vectorField> velocityCorrect
        (
            const vectorField&
        ) const override;

        //- Calculate the old boundary location
        virtual tmp<vectorField> oldBoundaryLocation() const override;

        virtual const pointField& oldPoints() const override
        {
            return oldPoints_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidBody&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
