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
    And the precice adapter solver object is based on the preCICE-
    adapter for OpenFOAM.

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
    (c) 2017-2023 Gerasimos Chourdakis
    (c) 2023 Engys Ltd.

Class
    Foam::fv::preciceAdapterSolverObject

Description
    This solver object prepares HELYX for preCICE coupling (precice.org).
    See the README.md for details.
    Look on precice.org for the latest version of the adapter.

    Example of solver object specification:
    \verbatim
    preciceAdapter
    {
        type preciceAdapter;
        libs ("libpreciceAdapter.so");
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: preciceAdapter | yes      |
    \endtable

SourceFiles
    preciceAdapterSolverObject.C

\*---------------------------------------------------------------------------*/

#ifndef preciceAdapterSolverObject_H
#define preciceAdapterSolverObject_H

#include "solverObjects/solverObject/solverObject.H"
#include "solverObjects/solverOption/SolverOption.H"

// OpenFOAM header files
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "clockValue/clockValue.H"

// Main adapter header file
#include "Adapter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class preciceAdapterSolverObject Declaration
\*---------------------------------------------------------------------------*/

class preciceAdapterSolverObject final
:
    public solverObject
{
    // Private data

// Some preCICE profiling
#ifdef ADAPTER_ENABLE_TIMINGS
    //- Cumulative time spent in the complete adapter function object (including preCICE calls)
    clockValue timeInAll_;
    clockValue timeInSetup_;
    clockValue timeInExecute_;
#endif

    //- Adapter's main class instance
    preciceAdapter::Adapter adapter_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("preciceAdapter");


    // Constructors

    //- Construct from components
    preciceAdapterSolverObject
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );


    //- No copy construct
    preciceAdapterSolverObject(preciceAdapterSolverObject&) = delete;

    //- No copy assignment
    void operator=(preciceAdapterSolverObject&) = delete;

    //- Destructor
    ~preciceAdapterSolverObject() final;


    // Solver Functions

        virtual void read(const dictionary& dict) override;

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual bool isFinalCorrector
        (
            const label corrector,
            const word& correctorName
        ) override;

        //- Return maximum acceptable time step
        virtual scalar getMaxTimeStep() override
        {
            return adapter_.getMaxTimeStep();
        }

        virtual void correct(const word&, const word&) override;

        //- Called at the bottom of the corrector loop
        virtual void endIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        ) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
