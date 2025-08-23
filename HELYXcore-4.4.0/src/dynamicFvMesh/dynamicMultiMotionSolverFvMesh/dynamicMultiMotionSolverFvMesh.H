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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::dynamicMultiMotionSolverFvMesh

Description
    Mesh motion described per cellZone. Individual motion solvers solve
    over whole domain but are only applied per cellZone.

SourceFiles
    dynamicMultiMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMultiMotionSolverFvMesh_H
#define dynamicMultiMotionSolverFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dynamicMultiMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMultiMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- The motion control function
        PtrList<motionSolver> motionPtr_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;


public:

    //- Runtime type information
    TypeName("dynamicMultiMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMultiMotionSolverFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicMultiMotionSolverFvMesh
        (
            const dynamicMultiMotionSolverFvMesh&
        ) = delete;


    //- Destructor
    ~dynamicMultiMotionSolverFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicMultiMotionSolverFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
