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
    (c) 2013-2024 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::displacementLaplacianPeriodicFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion displacement.

SourceFiles
    displacementLaplacianPeriodicFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLaplacianPeriodicFvMotionSolver_H
#define displacementLaplacianPeriodicFvMotionSolver_H

#include "fvMotionSolvers/displacement/laplacian/displacementLaplacianFvMotionSolver.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "helyxSolve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class helyxSolve;

/*---------------------------------------------------------------------------*\
             Class displacementLaplacianPeriodicFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLaplacianPeriodicFvMotionSolver
:
    public displacementLaplacianFvMotionSolver
{
    // Private data

        //- Amplification of first mode (t, amp)
        autoPtr<Function1<scalar>> amplification_;

        //- Time at max displacement
        scalar maxDisplacementTime_;

        //- Max displacement (for storing)
        mutable autoPtr<pointVectorField> maxPointDisplacement_;

        //- If maxPointDisplacement found at t=0, omit re-calc
        bool startFromDisp_;

        //- Create maxPointDisplacement field
        void createMaxDisplacementField();


public:

    //- Runtime type information
    TypeName("displacementLaplacianPeriodic");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementLaplacianPeriodicFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Construct from components
        displacementLaplacianPeriodicFvMotionSolver
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointVectorField& points0
        );

        //- Disallow default bitwise copy construction
        displacementLaplacianPeriodicFvMotionSolver
        (
            const displacementLaplacianPeriodicFvMotionSolver&
        ) = delete;


    //- Destructor
    ~displacementLaplacianPeriodicFvMotionSolver();


    // Member Functions

        //- Compute (max) displacement
        void computeCellDisplacement();
        void computePointDisplacement();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;
        tmp<pointField> curPoints(const scalar&) const;

        //- Solve for motion
        virtual void solve();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementLaplacianPeriodicFvMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
