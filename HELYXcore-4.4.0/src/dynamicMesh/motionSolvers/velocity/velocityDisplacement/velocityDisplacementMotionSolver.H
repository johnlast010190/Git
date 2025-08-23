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
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::velocityDisplacementMotionSolver

Description
    Mesh motion solver for a polyMesh. Wraps a displacement motion solver in a
    velocity motion solver.

SourceFiles
    velocityDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDisplacementMotionSolver_H
#define velocityDisplacementMotionSolver_H

#include "motionSolvers/velocity/velocityMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class displacementMotionSolver;

/*---------------------------------------------------------------------------*\
              Class velocityDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityDisplacementMotionSolver
:
    public velocityMotionSolver
{
    // Private data

        //- Displacement motion solver
        autoPtr<displacementMotionSolver> displacementMotionSolverPtr_;


    // Private Member Functions

        //- Get the boundary condition types for the point displacement
        wordList pointDisplacementBoundaryTypes() const;


public:

    //- Runtime type information
    TypeName("velocityDisplacement");


    // Constructors

        //- Construct from polyMesh and dictionary
        velocityDisplacementMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        velocityDisplacementMotionSolver
        (
            const velocityDisplacementMotionSolver&
        ) = delete;


    //- Destructor
    ~velocityDisplacementMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update geometry
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const velocityDisplacementMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
