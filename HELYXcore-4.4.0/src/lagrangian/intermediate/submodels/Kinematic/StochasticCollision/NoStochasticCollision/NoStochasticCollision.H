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
    Foam::NoStochasticCollision

Group
    grpLagrangianIntermediateCollisionSubModels

Description
    Dummy collision model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoStochasticCollision_H
#define NoStochasticCollision_H

#include "submodels/Kinematic/StochasticCollision/StochasticCollisionModel/StochasticCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NoStochasticCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoStochasticCollision
:
    public StochasticCollisionModel<CloudType>
{
protected:

    // Protected Member Functions

        //- Update the model
        virtual void collide(const scalar dt);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoStochasticCollision(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoStochasticCollision(const NoStochasticCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType>> clone() const
        {
            return autoPtr<StochasticCollisionModel<CloudType>>
            (
                new NoStochasticCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoStochasticCollision();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Kinematic/StochasticCollision/NoStochasticCollision/NoStochasticCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
