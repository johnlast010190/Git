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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2023 FOSS GP
    (c) 2024 Engys Ltd.

Class
    Foam::sinusFvMotionSolver

Description
    A demo motionSolver changing one spatial coordinate of all grid nodes based
    on the formula

    x = x0*(1 + amplitude*sin(2*pi/period + phase))

    x0: the initial mesh points coordinates
    amplitude: the amplitude of the change
    period: the period of the change
    phase: the phase of the change


SourceFiles
    sinusFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sinusMotionSolver_H
#define sinusMotionSolver_H

#include "motionSolvers/componentVelocity/componentVelocityMotionSolver.H"
#include "fvMotionSolvers/fvMotionSolver/fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sinusFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class sinusFvMotionSolver
:
    public componentVelocityMotionSolver,
    public fvMotionSolver
{
protected:

    // Protected data

        //- Period of the motion
        scalar period_;

        //- Radial frequency
        scalar radFrequency_;

        //- The motion phase
        scalar phase_;

        //- The amplitude multiplying the sin
        scalar amplitude_;

        //- Initial points
        pointField points0_;


private:

    // Private Member Functions

        //- No copy construct
        sinusFvMotionSolver
        (
            const sinusFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const sinusFvMotionSolver&) = delete;

public:

    //- Runtime type information
    TypeName("sinusFvMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        sinusFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~sinusFvMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        //  Does nothing
        virtual void solve(){};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
