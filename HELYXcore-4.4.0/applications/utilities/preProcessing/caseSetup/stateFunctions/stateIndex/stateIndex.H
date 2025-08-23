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
    (c) 2010-2024 Engys Ltd
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::stateIndex

Description
    Determines the type of solver based on dictionary settings and
    ouputs appropriate dictionaries.
    Also contains field prototypes and associated processing

SourceFiles
    stateIndex.C

\*---------------------------------------------------------------------------*/

#ifndef stateIndex_H
#define stateIndex_H

#include "containers/Lists/DynamicList/DynamicList.H"
#include "stateFunction/stateFunction.H"
#include "stateFunctions/global/globalState.H"
#include "regionModel/solutionRegions/solutionScheduler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class stateIndex Declaration
\*---------------------------------------------------------------------------*/

class stateIndex
{
    // Private data

        // Global references & dictionary input data

            //- Const reference to time
            const Time& runTime_;

            //- Region list
            HashTable<wordList> regionGroups_;

            solutionScheduler solutionSchedule_;

            //- Global input handling
            stateFunctions::globalState global_;


        // Regional input assembly

            //- Regional state data
            PtrList<stateFunction> regions_;


    // Private Member Functions

        //- Extract region information and output regionProperties file
        void regionProperties();

        //- Combine explicit, group and regExp input
        dictionary assembleRegionInput
        (
            const dictionary& input,
            const word& regionName,
            const word& rgID
        ) const;

        //- Compile input ditionaries
        void assembleDictionaries
        (
            const dictionary& input,
            dictionary& defaults
        );

        //- Load libraries from in-state controlDict
        void openLibs();


public:

    // Constructors

        //- Construct from dicitonary
        stateIndex
        (
            const Time& runTime,
            const dictionary& input,
            dictionary& defaults,
            const bool distributed,
            const bool collated
        );

        //- Disallow default bitwise copy construction
        stateIndex(const stateIndex&) = delete;


    // Member Functions

        //- Write dictionaries to file
        void writeDictionaries();

        //- Create mesh and store old fields
        void createRegionObjects();

        //- Modify mesh
        void modifyMesh(scalar writePause);

        //- Initialise fields
        void createFields();

        //- Switch reset boundaries
        void resetBoundaries(const bool& reset);

        const PtrList<stateFunction>& regions() const
        {
            return regions_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stateIndex&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
