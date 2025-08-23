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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::NoCollision

Group
    grpLagrangianIntermediateCollisionSubModels

Description
    Place holder for 'none' option

SourceFiles
    NoCollision.C

\*---------------------------------------------------------------------------*/

#ifndef NoCollision_H
#define NoCollision_H

#include "submodels/Kinematic/CollisionModel/CollisionModel/CollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class NoCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoCollision
:
    public CollisionModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        NoCollision(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoCollision(const NoCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType>> clone() const
        {
            return autoPtr<CollisionModel<CloudType>>
            (
                new NoCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoCollision();


    // Member Functions

        //- Return the number of times to subcycle the current
        //  timestep to meet the criteria of the collision model.  For
        //  this model this will always be 1.
        virtual label nSubCycles() const;

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Indicates whether model determines wall collisions or not,
        //  used to determine what value to use for wallImpactDistance
        virtual bool controlsWallInteraction() const;

        // Collision function
        virtual void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Kinematic/CollisionModel/NoCollision/NoCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
