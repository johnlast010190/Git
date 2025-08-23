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

#ifndef solutionInstance_H
#define solutionInstance_H

#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/lists/hashedWordList.H"
#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionInstance
:
    public solutionInstanceRegistry
{
protected:

    PtrList<fvSolutionRegistry> solutionRegistryPtrList_;
    UPtrList<fvMesh> meshPtrList_;
    bool backwardCompatibility_;

public:


    // Constructors

        //- Construct from components
        solutionInstance
        (
            const Time& runTime,
            const List<List<word>>& instanceRegions,
            const dictionary& dict,
            const word& name,
            bool backwardCompatibility = false
        );

        //- Disallow default bitwise copy construction
        solutionInstance(const solutionInstance&) = delete;


    //- Destructor
    ~solutionInstance() override;

    // Member Functions

        inline const dictionary& dict() const
        {
            return  solutionMeshDict_;
        }

        const List<wordList>& solutionRegions() const
        {
            return solutionRegions_;
        }

        const fvSolutionRegistry& getSolutionRegistry(const Time& runTime);

        const PtrList<fvSolutionRegistry>& getSolutionRegistryList
        (
            const Time& runTime
        );

        const PtrList<fvSolutionRegistry>& getSolutionRegistryListOldFormat
        (
            const Time& runTime
        );

        fvMesh* lookupOrCreateMesh
        (
            const Time& runTime,
            const word& meshName
        );

        void releaseInstaceFields();

        void releaseUnusedMeshes(const Time& runTime, const word& meshName);

        void deactivate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solutionInstance&) = delete;
};

}

#endif
