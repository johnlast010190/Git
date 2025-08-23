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
    (c) 2024 Engys Ltd.

Class
    Foam::fv::waves

Group
    grpFvOptionsSources

Description
    Wave relaxation zones

SourceFiles
    waves.C

\*---------------------------------------------------------------------------*/

#ifndef waves_H
#define waves_H

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
              Class waves Declaration
\*---------------------------------------------------------------------------*/

class waves
:
    public solverObject
{
protected:

    //- Multiphase thermo model
    multiphaseThermo* thermoPtr_;

    //- The first non-passive volume fraction field
    volScalarField* alpha1_;

    //- waves2Foam relaxation zones
    autoPtr<relaxationZone> relaxing_;

public:

    //- Runtime type information
    TypeName("waves");


    // Constructors

        //- Construct from components
        waves
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        waves(waves&) = delete;

    void operator=(waves&) = delete;


    //- Destructor
    virtual ~waves()
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
