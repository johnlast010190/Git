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
    Foam::deformingBodyMotionFunctions::staticMotion

Description
    Motion for burning interface.
    This class calculates the interface velocity and velocity boudnary values

SourceFiles
    staticMotion.C

\*---------------------------------------------------------------------------*/

#ifndef staticMotion_H
#define staticMotion_H

#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "fields/Fields/vectorField/vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace deformingBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class staticMotion Declaration
\*---------------------------------------------------------------------------*/

class staticMotion
:
    public deformingBodyMotionFunction
{
    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> normalVelocity(const label&) const;


public:

    //- Runtime type information
    TypeName("staticMotion");


    // Constructors

        //- Construct from components
        staticMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        staticMotion(const staticMotion&) = delete;


    //- Destructor
    virtual ~staticMotion();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary&);

        //- update
        virtual void update();

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const staticMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace deformingBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
