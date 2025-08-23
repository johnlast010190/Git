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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::motionFunction

Description
    Base class for defining interface motion for GIB.

SourceFiles
    deformingBodyMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef deformingBodyMotionFunction_H
#define deformingBodyMotionFunction_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
              Class fvMeshGIBChangers::motionFunction Declaration
\*---------------------------------------------------------------------------*/

class motionFunction
{
protected:

    // Protected Data

        const fvMesh& mesh_;

        dictionary DBMFCoeffs_;

        const Time& time_;


public:

    //- Runtime type information
    TypeName("motionFunction");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        motionFunction,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& SBMFCoeffs,
            const Time& runTime
        ),
        (mesh, SBMFCoeffs, runTime)
    );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        motionFunction
        (
            const fvMesh&,
            const dictionary&,
            const Time&
        );

        //- Disallow default bitwise copy construction
        motionFunction(const motionFunction&) = delete;


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<motionFunction> New
        (
            const fvMesh&,
            const dictionary&,
            const Time&
        );


    //- Destructor
    virtual ~motionFunction();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;

        //- Update
        virtual void update() = 0;

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfacePointsVelocity(const label&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const motionFunction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
