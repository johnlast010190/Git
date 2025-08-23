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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef solutionScheduler_H
#define solutionScheduler_H

#include "db/Time/Time.H"
#include "global/argList/argList.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/lists/hashedWordList.H"
#include "regionModel/solutionInstance/solutionInstance.H"
#include "regionModel/solutionRegionMapping/solutionRegionMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionScheduler
{

    const Time& runTime_;


protected:

    IOdictionary dict_;

    solutionRegionMapping regionMap_;

    wordList instanceNames_;

    hashedWordList meshNames_;

    HashSet<word> selectedInstances_;

    label numberOfInstances_;

    label instanceCounter_;

    label currentInstanceIndex_;

    label nextInstanceIndex_;

    PtrList<solutionInstance> solutionInstancePtrList_;

    bool defaultBehavior_;

    dictionary solutionScheduleDict_;

    dictionary solutionMeshDict_;

    void assignUserInputValues();

    void selectOnly(const HashSet<word>& selectedInstances);

    void checkSolutionMesh();

    void buildInstances();

    void checkValidInput();

    void assignDefaultValues();


public:

    // Constructors

        //- Construct from time
        explicit solutionScheduler
        (
            argList& args,
            const Time& runTime
        );

        //- Construct from time, dict and flag
        solutionScheduler
        (
            const Time& runTime,
            const IOdictionary& IOdict,
            bool caseSetup = false
        );

        //- Disallow default bitwise copy construction
        solutionScheduler(const solutionScheduler&) = delete;


    //- Destructor
    ~solutionScheduler();


    // Member Functions

        bool loop();

        const fvSolutionRegistry& getActiveSolutionRegistry(const Time& runTime);

        const PtrList<fvSolutionRegistry>& getActiveSolutionRegistryList
        (
            const Time& runTime
        );

        dictionary getSchedulerDictionary() const;

        dictionary getSolutionMeshDictionary() const;

        inline bool active() const
        {
            return !defaultBehavior_;
        }

        const PtrList<solutionInstance>& instances() const
        {
            return solutionInstancePtrList_;
        }

        const solutionInstance& getCurrentInstance() const
        {
            return solutionInstancePtrList_[currentInstanceIndex_];
        }


    // Member Operators

        void operator++();

        void operator++(int);

        //- Disallow default bitwise assignment
        void operator=(const solutionScheduler&) = delete;
};

}

#endif
