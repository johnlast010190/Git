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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::fv::turbulenceSolver

Description
    Corrects the turbulence model. Creates compressible model on the object
    registry if none found.

SourceFiles
    turbulenceSolver.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceSolver_H
#define turbulenceSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "turbulenceModel.H"
#include "transportModel/transportModel.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianMultiphaseSystem;

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class turbulenceSolver Declaration
\*---------------------------------------------------------------------------*/

class turbulenceSolver
:
    public solverObject
{
    // Private data

        //- The turbulence model pointer, either compressible::turbulenceModel
        //  or incompressible::turbulenceModel
        turbulenceModel* turbulencePtr_;

        //- Pointer to phase system for multiple phasic turbulence models
        eulerianMultiphaseSystem* phaseSystemPtr_;

        //- Laminar transport model if incompressible
        transportModel* laminarTransportPtr_;

        //- Solution control
        solutionControl* solnControlPtr_;

        //- single or multi-region
        bool singleRegion_;

public:

    //- Runtime type information
    TypeName("turbulenceSolver");


    // Constructors

        //- Construct from components
        turbulenceSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        turbulenceSolver(turbulenceSolver&) = delete;

    void operator=(turbulenceSolver&) = delete;


    // Solver Functions

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word& solveName, const word&) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
