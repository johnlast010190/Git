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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::dynamicMotionSolverFvMesh

Description
    The dynamicMotionSolverFvMesh

SourceFiles
    dynamicMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverFvMesh_H
#define dynamicMotionSolverFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                           Class dynamicMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        autoPtr<motionSolver> motionPtr_;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicMotionSolverFvMesh(const dynamicMotionSolverFvMesh&) = delete;


    //- Destructor
    ~dynamicMotionSolverFvMesh();


    // Member Functions

        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Write the underlying polyMesh and other data
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
