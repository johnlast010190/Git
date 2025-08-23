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
    (c) 2019-2022 Engys Ltd.

Class
    Foam::coordinateFrame

Description
    Base class for reference frames. Handles information on position,
    velocity and acceleration of the frame and provides to sub-systems.

\*---------------------------------------------------------------------------*/

#ifndef coordinateFrame_H
#define coordinateFrame_H

#include "meshes/MeshObject/MeshObject.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fvMesh/fvMesh.H"
#include "coordinate/systems/coordinateSystem.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "referenceFrames/coordinateFrameRegistry/coordinateFrameRegistry.H"
#include "referenceFrames/coordinateFrameState/coordinateFrameState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coordinateFrame Declaration
\*---------------------------------------------------------------------------*/

class coordinateFrame
:
    public MeshObject<fvMesh, UpdateableMeshObject, coordinateFrame>,
    public coordinateFrameState
{

    // Private Data

        //- Pointer list to all parents
        mutable UPtrList<coordinateFrame> parents_;

        //- Global septerninon transformation
        mutable septernion globalTransform_;


    // Private member functions

        //- Initialise pointers to all parents
        void initParents() const;


protected:

        mutable coordinateFrameRegistry coorFrameReg_;

        //- Coordinate Frame Dictionary
        dictionary frameDict_;

        //- Parent coordinate frame
        mutable coordinateFrame* parentFramePtr_;

        //- Name of parent frame
        word parentFrameName_;


public:

    //- Runtime type information
    TypeName("coordinateFrame");


    // Constructors

        //- Construct from components
        coordinateFrame
        (
            const fvMesh& mesh,
            const dictionary& extRefDict,
            const word& frameName
        );


    //- Destructor
    virtual ~coordinateFrame();

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateFrame,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& frameName
        ),
        (mesh, dict, frameName)
    );


    // Selectors

        //- Constructor to create new frame or return reference to the old one
        static coordinateFrame& New
        (
            const fvMesh& mesh,
            const word& frameName
        );

        //- Constructor to create new frame or return reference to the old one
        static coordinateFrame* New
        (
            const fvMesh& mesh,
            const word& frameName,
            const dictionary& dict
        );

        //- Helper constructor to omit lookup "referenceFrame"
        static coordinateFrame* lookupNew
        (
            const fvMesh& mesh,
            const dictionary& dictWithFrameName
        );

        //- Construct and return global frame
        static coordinateFrame* globalFrame(const fvMesh& mesh);


    // Static Functions

        //- Update stat of all reference frames in this object registry
        static void updateStates(const objectRegistry& obr)
        {
            const wordList frameNames(obr.names<coordinateFrame>());
            forAll(frameNames, namei)
            {
                obr.lookupObjectRef<coordinateFrame>
                (
                    frameNames[namei]
                ).updateState();
            }
        }

        //- Update stat of all reference frames in this object registry
        static void resetUpdates
        (
            const objectRegistry& obr,
            label outerCorrectorUpdate = 1
        )
        {
            const wordList frameNames(obr.names<coordinateFrame>());
            forAll(frameNames, namei)
            {
                obr.lookupObjectRef<coordinateFrame>
                (
                    frameNames[namei]
                ).resetUpdate(outerCorrectorUpdate);
            }
        }


    // Mesh motion and mophing

        //- Reference to coordinate frame registry
        const coordinateFrameRegistry& registry()
        {
            return coorFrameReg_;
        };

        //- Is mesh moving
        bool moving() const
        {
            return mesh().moving();
        }

        //- Update after mesh motion
        bool movePoints() override
        {
            return true;
        }

        //- Return reference to parents
        UPtrList<coordinateFrame>& parents() const
        {
            initParents();
            return parents_;
        }

        //- Any of the frame and its parents is moving
        bool anyDynamic() const
        {
            forAll(parents(), i)
            {
                if (parents()[i].isDynamic())
                {
                    return true;
                }
            }
            return false;
        }

        //- Any of the frame and its parents is moving
        bool anyIncremental() const
        {
            forAll(parents(), i)
            {
                if (parents()[i].isIncrementalMotion())
                {
                    return true;
                }
            }
            return false;
        }

        //- Update after topo change
        void topoChange(const polyTopoChangeMap&) override
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) override
        {}

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&) override
        {}


    // Access Functions

        //- Return parent frame name
        const word& parentFrameName() const
        {
            return parentFrameName_;
        }

        //- Parent frame access
        const coordinateFrame& parentFrame() const
        {
            return *parentFramePtr_;
        }

        //- Is parent frame valid
        bool validParentFrame() const;

        //- Add parent frame (only if it is valid)
        bool addParentFrame(bool addParentFrame) const
        {
            return validParentFrame() && addParentFrame;
        }

        //- Current rotational speed
        virtual vector Omega() const
        {
            updateState();
            return velocity().second();
        }

        //- Center of rotation
        virtual const vector& CofR() const
        {
            return coorSys().origin();
        }

        //- Initial center of rotation
        virtual const vector& CofR0() const
        {
            return coorSys0().origin();
        }

        //- Return coordinate system in local definition
        autoPtr<coordinateSystem> localCoordSystem() const;

        //- Update state of reference frame
        virtual void updateState(bool construction = false) const;

        //- Return frame relative velocity vector from position(s)
        virtual tmp<volVectorField> frameVelocity
        (
            const volVectorField& positions,
            bool addParentFrames = true
        ) const;

        //- Return velocity vector field from position(s)
        virtual tmp<vectorField> frameVelocity
        (
            const vectorField& positions,
            bool addParentFrames = true
        ) const;

        //- Return frame velocity for specific point
        virtual vector frameVelocity
        (
            const vector& position,
            bool addParentFrames = true
        ) const;

        //- Calculating velocity contributions from reference frames
        //  to the boundary patch. It should handle MRF and dynamic meshes
        //  working together.
        virtual void calculateBoundaryVelocity
        (
            vectorField& Up,
            const label patchi,
            const objectRegistry& obr,
            const word& UName,
            bool inletFlux
        ) const;

        //- Add all normal components of velocity, however replace mesh motion
        //  with normal recomputed through the mesh flux
        virtual void inletFluxVelocity
        (
            vectorField& Up,
            const label patchi,
            const objectRegistry& obr,
            const word& UName,
            bool addedMeshComponent = false
        ) const;

        //- Attaching patch to the MRF
        virtual void attachPatch(const label& patchi) const;

        //- Check if the patch is attached to the MRF
        virtual bool isAttachToMRF(const label& patchi) const;

        //- Return the frame motion transformation septernion
        virtual const septernion& transformation
        (
            bool addParentFrames = true,
            label nCorr = -1
        ) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif


// ************************************************************************* //
