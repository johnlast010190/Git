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
    (c) 2025 Engys Ltd.

Class
    Foam::fv::phaseChangeSource

Group
    grpFvOptionsSources

Description
    Source terms for wave relaxation zones

SourceFiles
    phaseChangeSource.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeSource_H
#define phaseChangeSource_H

#include "solverObjects/solverObject/solverObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class relaxationZone;
class multiphaseThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
              Class phaseChangeSource Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeSource
:
    public solverObject
{
protected:

    //- Multiphase thermo model
    multiphaseThermo* thermoPtr_;

public:

    //- Runtime type information
    TypeName("phaseChangeSource");


    // Constructors

        //- Construct from components
        phaseChangeSource
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        phaseChangeSource(phaseChangeSource&) = delete;

    void operator=(phaseChangeSource&) = delete;


    //- Destructor
    virtual ~phaseChangeSource()
    {}


    // Solver Functions

        virtual bool initialise() override;

        virtual void read(const dictionary& dict) override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;

    // Source term functions

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
        ) override;

        //- Add contribution to p and alpha equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Add contribution to p equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
