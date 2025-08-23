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
    Foam::ORourkeCollision

Description
    Collision model by P.J. O'Rourke.


\*---------------------------------------------------------------------------*/

#ifndef ORourkeCollision_H
#define ORourkeCollision_H

#include "submodels/Kinematic/StochasticCollision/StochasticCollisionModel/StochasticCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class ORourkeCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ORourkeCollision
:
    public StochasticCollisionModel<CloudType>
{
protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        const liquidMixtureProperties& liquids_;

        //- Coalescence activation switch
        Switch coalescence_;


    // Protected Member Functions

        //- Main collision routine
        virtual void collide(const scalar dt);

        //- Collide parcels and return true if mass has changed
        virtual bool collideParcels
        (
            const scalar dt,
            parcelType& p1,
            parcelType& p2,
            scalar& m1,
            scalar& m2
        );

        // 1 is the larger drop and 2 is the smaller
        virtual bool collideSorted
        (
            const scalar dt,
            parcelType& p1,
            parcelType& p2,
            scalar& m1,
            scalar& m2
        );


public:

    //- Runtime type information
    TypeName("ORourke");


    // Constructors

        //- Construct from dictionary
        ORourkeCollision
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& modelName = typeName
        );

        //- Construct copy
        ORourkeCollision(const ORourkeCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType>> clone() const
        {
            return autoPtr<StochasticCollisionModel<CloudType>>
            (
                new ORourkeCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ORourkeCollision();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/StochasticCollision/ORourkeCollision/ORourkeCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
