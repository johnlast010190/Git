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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::cellSet

Description
    A collection of cell labels.

SourceFiles
    cellSet.C

\*---------------------------------------------------------------------------*/

#ifndef cellSet_H
#define cellSet_H

#include "sets/topoSets/topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellSet Declaration
\*---------------------------------------------------------------------------*/

class cellSet
:
    public topoSet
{
public:

    //- Runtime type information
    TypeName("cellSet");


    // Constructors

        //- Construct from IOobject. No checking.
        cellSet(const IOobject& obj);

        //- Construct from polyMesh and name. Checks for valid cell ids.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct empty from size of labelHashSet
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const label size,
            writeOption w=NO_WRITE
        );

        //- Construct from existing set
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            writeOption w=NO_WRITE
        );

        //- Construct from labelHashSet
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet& set,
            writeOption w=NO_WRITE
        );

        //- Construct from additional list of labels for the labelHashSet
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const UList<label>& set,
            writeOption w=NO_WRITE
        );

        //- Disallow default bitwise copy construction
        cellSet(const cellSet&) = delete;


        // Used for tetMesh cellSet only.

            //- Construct from objectRegistry and name.
            cellSet
            (
                const Time&,
                const word& name,
                readOption r=MUST_READ,
                writeOption w=NO_WRITE
            );

            //- Construct empty from objectRegistry.
            cellSet
            (
                const Time&,
                const word& name,
                const label size,
                writeOption w=NO_WRITE
            );

            //- Construct from labelHashSet
            cellSet
            (
                const Time&,
                const word& name,
                const labelHashSet&,
                writeOption w=NO_WRITE
            );


    //- Destructor
    virtual ~cellSet();


    // Member functions

        //- Sync cellSet across coupled patches.
        virtual void sync(const polyMesh& mesh)
        {}

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels.
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const polyDistributionMap&);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
