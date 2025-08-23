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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::adjointMotion

Description
    Deforming body adjoint motion function.

SourceFiles
    adjointMotion.C

\*---------------------------------------------------------------------------*/

#ifndef adjointMotion_H
#define adjointMotion_H

#include "deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
               Class fvMeshGIBChangers::adjointMotion Declaration
\*---------------------------------------------------------------------------*/

class adjointMotion
:
    public motionFunction
{
    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> normalVelocity(const label&) const;


public:

    //- Runtime type information
    TypeName("adjointMotion");


    // Constructors

        //- Construct from components
        adjointMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~adjointMotion();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary& DBMFCoeffs) override
        {
            return true;
        }

        //- Update
        virtual void update() override
        {}

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const override;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adjointMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
