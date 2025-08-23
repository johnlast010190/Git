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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::sixDoFRigidBodyMotionRestraints::tabulatedAxialAngularSpring

Description
    sixDoFRigidBodyMotionRestraints model.  Axial angular spring with moment
    values drawn from an interpolation table.  Linear damping.

SourceFiles
    tabulatedAxialAngularSpring.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedAxialAngularSpring_H
#define tabulatedAxialAngularSpring_H

#include "sixDoFRigidBodyMotion/restraints/sixDoFRigidBodyMotionRestraint/sixDoFRigidBodyMotionRestraint.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "interpolations/interpolationTable/interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedAxialAngularSpring Declaration
\*---------------------------------------------------------------------------*/

class tabulatedAxialAngularSpring
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Reference orientation where there is no moment
        tensor refQ_;

        //- Global unit axis around which the motion is sprung
        vector axis_;

        //- Spring moment interpolation table, depending on angleFormat
        interpolationTable<scalar> moment_;

        //- Boolean stating whether the angle around the axis needs to
        //  be converted to degrees before asking the
        //  interpolationTable for a moment value
        bool convertToDegrees_;

        //- Damping coefficient (Nms/rad)
        scalar damping_;


public:

    //- Runtime type information
    TypeName("tabulatedAxialAngularSpring");


    // Constructors

        //- Construct from components
        tabulatedAxialAngularSpring
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new tabulatedAxialAngularSpring(*this)
            );
        }


    //- Destructor
    virtual ~tabulatedAxialAngularSpring();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFRigidBodyMotionRestraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
