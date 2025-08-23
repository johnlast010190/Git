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
    (c) 2018-2024 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::burningMotion

Description
    Motion for burning interface.

    This class calculates the interface velocity and velocity boundary values.

SourceFiles
    burningMotion.C

\*---------------------------------------------------------------------------*/

#ifndef burningMotion_H
#define burningMotion_H

#include "deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunction.H"
#ifndef HELYX_SWAK4FOAM_DISABLED
#include "expressionSource.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
               Class fvMeshGIBChangers::burningMotion Declaration
\*---------------------------------------------------------------------------*/

class burningMotion
:
    public motionFunction
{
    // Private Data

        //- Motion type
        word mode_;

        //- Solid density
        scalar rhos_;

        //- Constants
        autoPtr<Function1<scalar>> uniformBurnRate_;

        //- User-defined field name
        word fieldName_;

        //- Expression entered by user
        #ifndef HELYX_SWAK4FOAM_DISABLED
        mutable autoPtr<expressionSource<scalar>> expression_;
        #endif

        //- Scale factor to accelerate the burn progression without
        //  influencing boundary velocity
        scalar scaleFactor_;


    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> burnRate(const label&) const;


public:

    //- Runtime type information
    TypeName("burningMotion");


    // Constructors

        //- Construct from components
        burningMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        burningMotion(const burningMotion&) = delete;


    //- Destructor
    virtual ~burningMotion();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary&) override;

        //- Update
        virtual void update() override
        {}

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const override;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const burningMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
