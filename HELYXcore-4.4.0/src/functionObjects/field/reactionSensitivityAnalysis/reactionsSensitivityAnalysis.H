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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::reactionsSensitivityAnalysis

Group
    grpFieldFunctionObjects grpThermophysicalFunctionObjects

Description
    This function object creates four data files named:

    - "consumption"    :   consumption rate
    - "production"     :   destruction rate
    - "productionInt"  :   integral between dumps of the production rate
    - "consumptionInt" :   integral between dumps of the consumption rate

    The function object indicates reaction rates of creation or destruction
    of species in each reaction.


SourceFiles
    reactionsSensitivityAnalysis.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjescts_reactionsSensitivityAnalysis_H
#define functionObjescts_reactionsSensitivityAnalysis_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFieldsFwd.H"
#include "basicChemistryModel/basicChemistryModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class reactionsSensitivityAnalysis Declaration
\*---------------------------------------------------------------------------*/

template<class chemistryType>
class reactionsSensitivityAnalysis
:
    public fvMeshFunctionObject,
    public writeFile

{
    // Private data

        //- List list for species production
        scalarListList production_;

        //- List list for species consumption
        scalarListList consumption_;

        //- List list for species production integral
        scalarListList productionInt_;

        //- List list for species consumption integral
        scalarListList consumptionInt_;

        //- Start time of integration
        scalar startTime_;

        //- End time of integration
        scalar endTime_;

        //- Word list of species
        wordList speciesNames_;

        //-Number of reactions
        label nReactions_;


         // File streams

            //- Integrated coefficients
            autoPtr<OFstream> prodFilePtr_;

            //- Moment coefficient
            autoPtr<OFstream> consFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> prodIntFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> consIntFilePtr_;



    // Private Member Functions


        //- Create file names for forces and bins
        void createFileNames();

        //- Output file header information
        void writeFileHeader(OFstream& os);

        //- Calculate production and destruction of each species
        void calculateSpeciesRR(const basicChemistryModel&);

        //- Write species production/consumption rates
        void writeSpeciesRR();


public:

    //- Runtime type information
    TypeName("reactionsSensitivityAnalysis");


    // Constructors

        //- Construct from Time and dictionary
        reactionsSensitivityAnalysis
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        reactionsSensitivityAnalysis
        (
            const reactionsSensitivityAnalysis&
        ) = delete;


    //- Destructor
    virtual ~reactionsSensitivityAnalysis();


    // Member Functions

        //- Read the reactionsSensitivityAnalysis data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the reactionsSensitivityAnalysis and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactionsSensitivityAnalysis&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "reactionSensitivityAnalysis/reactionsSensitivityAnalysis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
