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
    (c) 2016 OpenCFD Ltd
    (c) 2024 Engys Ltd.

Class
    Foam::dynamicMotionSolverTopoFvMesh

Description
    Sample topoChangerFvMesh that moves an object using a motion solver

SourceFiles
    dynamicMotionSolverTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverTopoFvMesh_H
#define dynamicMotionSolverTopoFvMesh_H

#include "topoChangerFvMesh/topoChangerFvMesh.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class dynamicMotionSolverTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Motion solver
        autoPtr<motionSolver> motionPtr_;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit dynamicMotionSolverTopoFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicMotionSolverTopoFvMesh
        (
            const dynamicMotionSolverTopoFvMesh&
        ) = delete;


    //- Destructor
    virtual ~dynamicMotionSolverTopoFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverTopoFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
