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
    (c) 2024 Engys Ltd.

Class
    Foam::multiSolidBodyMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    multiSolidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBodyMotionSolver_H
#define multiSolidBodyMotionSolver_H

#include "motionSolvers/displacement/points/pointsMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSolidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyMotionSolver
:
    public pointsMotionSolver
{
    // Private data

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;


public:

    //- Runtime type information
    TypeName("multiSolidBodyMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBodyMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        multiSolidBodyMotionSolver(const multiSolidBodyMotionSolver&) = delete;


    //- Destructor
    ~multiSolidBodyMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const override;

        //- Solve for motion
        virtual void solve() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSolidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
