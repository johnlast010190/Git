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

Class
    Foam::RBD::joints::Rs

Group
    grpRigidBodyDynamicsJoints

Description
    Spherical joint for rotation about the x/y/z-axes using a quaternion
    (Euler parameters) to avoid gimble-lock.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
        Chapter 4.
    \endverbatim

SourceFiles
    Rs.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_joints_Rs_H
#define RBD_joints_Rs_H

#include "joints/joint/joint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{

/*---------------------------------------------------------------------------*\
                         Class Rs Declaration
\*---------------------------------------------------------------------------*/

class Rs
:
    public joint
{

public:

    //- Runtime type information
    TypeName("Rs");


    // Constructors

        //- Construct for given model
        Rs();

        //- Construct for given model from dictionary
        Rs(const dictionary& dict);

        //- Clone this joint
        virtual autoPtr<joint> clone() const;


    //- Destructor
    virtual ~Rs();


    // Member Functions

        //- Return true as this joint describes rotation using a quaternion
        virtual bool unitQuaternion() const;

        //- Update the model state for this joint
        virtual void jcalc
        (
            joint::XSvc& J,
            const scalarField& q,
            const scalarField& qDot
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace joints
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
