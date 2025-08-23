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

#ifndef solutionInstanceRegistry_H
#define solutionInstanceRegistry_H

#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionInstanceRegistry
:
    public objectRegistry
{
protected:

    dictionary solutionMeshDict_;

    DynamicList<word> meshNames_;

    List<wordList> solutionRegions_;

    wordList regionNames_;

    mutable bool active_;

    //- Flag to tell which instances are selected for the scheduler
    mutable bool selected_;

    label getTotalNumberOfRegions();

    void generateMeshList();

    label whichRegionIndex(const word& regionName) const;


public:

    // Constructors

        //- Construct from components
        solutionInstanceRegistry
        (
            const objectRegistry& obr,
            const List<List<word>>& instanceRegions,
            const dictionary& dict,
            const word& name
        );

        //- Disallow default bitwise copy construction
        solutionInstanceRegistry(const solutionInstanceRegistry&) = delete;

    //- Destructor
    ~solutionInstanceRegistry();


    // Member Functions

        const word whichMesh(const word& regionName) const;

        inline const bool& isActive() const
        {
            return active_;
        }

        inline const bool& isSelected() const
        {
            return selected_;
        }

        inline const wordList& regionNames() const
        {
            return regionNames_;
        }

        inline const wordList& meshNames() const
        {
            return meshNames_;
        }

        void unselectInstance()
        {
            selected_ = false;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solutionInstanceRegistry&) = delete;
};

}

#endif
