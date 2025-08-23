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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::RBD::rigidBodyMotion

Description
    Six degree of freedom motion for a rigid body.

    Angular momentum stored in body fixed reference frame.  Reference
    orientation of the body (where Q = I) must align with the cartesian axes
    such that the Inertia tensor is in principle component form.  Can add
    restraints (e.g. a spring) and constraints (e.g. motion may only be on a
    plane).

    The time-integrator for the motion is run-time selectable with options for
    symplectic (explicit), Crank-Nicolson and Newmark schemes.

SourceFiles
    rigidBodyMotionI.H
    rigidBodyMotion.C
    rigidBodyMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyMotion_H
#define rigidBodyMotion_H

#include "rigidBodyModel/rigidBodyModel.H"
#include "rigidBodyModelState/rigidBodyModelState.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

// Forward declarations
class rigidBodySolver;

/*---------------------------------------------------------------------------*\
                      Class rigidBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyMotion
:
    public rigidBodyModel
{
    friend class rigidBodySolver;

    // Private data

        //- Motion state data object
        rigidBodyModelState motionState_;

        //- Motion state data object for previous time-step
        rigidBodyModelState motionState0_;

        //- Initial transform for external forces to the bodies reference frame
        List<spatialTransform> X00_;

        //- Acceleration relaxation coefficient
        scalar aRelax_;

        //- Acceleration damping coefficient (for steady-state simulations)
        scalar aDamp_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;

        //- Motion solver
        autoPtr<rigidBodySolver> solver_;


    // Private Member Functions

        //- Initialize the body-state
        void initialize();


public:

    // Constructors

        //- Construct null
        rigidBodyMotion();

        //- Construct from dictionary
        rigidBodyMotion
        (
            const dictionary& dict
        );

        //- Construct from constant and state dictionaries
        rigidBodyMotion
        (
            const dictionary& dict,
            const dictionary& stateDict
        );

        //- Disallow default bitwise copy construction
        rigidBodyMotion(const rigidBodyMotion&) = delete;


    //- Destructor
    ~rigidBodyMotion();


    // Member Functions

        // Access

            //- Return the report Switch
            inline bool report() const;

            //- Return the motion state
            inline const rigidBodyModelState& state() const;

            //- Return the motion state for modification
            inline rigidBodyModelState& state();

            //- Return the initial transform to the global frame for the
            //  given body
            spatialTransform X00(const label bodyId) const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void newTime();


        // Update state

            //- Calculate and optionally relax the joint acceleration qDdot from
            //  the joint state q, velocity qDot, internal force tau (in the
            //  joint frame) and external force fx (in the global frame)
            void forwardDynamics
            (
                rigidBodyModelState& state,
                const scalarField& tau,
                const Field<spatialVector>& fx
            ) const;

            //- Integrate velocities, orientation and position
            //  for the given time-step
            void solve
            (
                scalar deltaT,
                const scalarField& tau,
                const Field<spatialVector>& fx
            );

            //- Report the status of the motion of the given body
            void status(const label bodyID) const;


        // Transformations

            //- Return the transformation of bodyID relative to the initial time
            spatialTransform transform0(const label bodyID) const;


        //- Write
        void write(Ostream&) const;

        //- Read coefficients dictionary and update system parameters,
        //  constraints and restraints but not the current state
        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBodyMotion/rigidBodyMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
