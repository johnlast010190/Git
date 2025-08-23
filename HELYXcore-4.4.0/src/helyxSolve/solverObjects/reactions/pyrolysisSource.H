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
    Foam::fv::pyrolysisSource

Group
    grpFvOptionsSources

Description
    Thermal and species source terms from combustion model.

SourceFiles
    pyrolysisSource.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisSource_H
#define pyrolysisSource_H

#include "solverObjects/solverObject/solverObject.H"
#include "combustionModel/combustionModel.H"
#include "pyrolysisModels/pyrolysisModel/pyrolysisModelCollection.H"
#include "solidChemistryModel/solidChemistryModel.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class pyrolysisSource Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisSource
:
    public solverObject
{
protected:


    // Protected Member Data

        //- Reaction model
        autoPtr<regionModels::pyrolysisModelCollection> pyrolysisPtr_;

        //- Pointer to solution control
        solutionControl* solnControlPtr_;

        //- Solid region diffusion number
        scalar diNum_;

        //- Maximum diffusion number
        scalar maxDi_;


public:

    //- Runtime type information
    TypeName("pyrolysisSource");


    // Constructors

        //- Construct from components
        pyrolysisSource
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        pyrolysisSource(pyrolysisSource&) = delete;


    //- Destructor
    virtual ~pyrolysisSource()
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

        scalar getMaxTimeStep() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(pyrolysisSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
