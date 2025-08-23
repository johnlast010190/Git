/*---------------------------------------------------------------------------* \
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
    (c) 2019-2022 Engys Ltd.
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fv::fixedTemperature

Group
    grpFvOptionsSolvers

Description
    Sets the fluid energy such as to ensure the temperature distribution does
    not change.

SourceFiles
    fixedTemperature.C

SeeAlso
    Foam::fv::solidHeatConductionSolver

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperature_H
#define fixedTemperature_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class fixedTemperature Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperature
:
    public solverObject
{

protected:

    basicThermo* thermoPtr_;

public:

//- Runtime type information
    TypeName("fixedTemperature");

// Constructors

    // - Construct from components
    fixedTemperature
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

    fixedTemperature(fixedTemperature&) = delete;

    void operator=(fixedTemperature&) = delete;

//- Destructor
    virtual ~fixedTemperature()
    {}

// Member Functions

    // Solver functions

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& derivedFields,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
