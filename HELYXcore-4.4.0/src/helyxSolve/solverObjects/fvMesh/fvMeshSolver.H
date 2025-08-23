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
    (c) 2021-2024 Engys Ltd

Class
    Foam::fv::fvMeshSolver

Description
    Solver object for updating the fvMesh

SourceFiles
    fvMeshSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshSolver_H
#define fvMeshSolver_H

#include "solverObjects/solverObject/solverObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSolver
:
    public solverObject
{
public:

    //- Runtime type information
    TypeName("fvMeshSolver");

    // Constructors

        //- Construct from components
        fvMeshSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Dissallow bitwise copy construction
        fvMeshSolver(fvMeshSolver&) = delete;


    // Member Functions

        //- Initialise the solver object
        virtual bool initialise() override;

        //- Get the solve graph
        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        //- Return maximum acceptable time step
        virtual scalar getMaxTimeStep() override;

        //- Correct
        virtual void correct(const word&, const word&) override;


    // Member Operators

        //- Dissallow bitwise assignment
        void operator=(fvMeshSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
