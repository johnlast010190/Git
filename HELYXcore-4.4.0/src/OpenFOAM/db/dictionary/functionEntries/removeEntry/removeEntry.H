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
    (c) 2024 Engys Ltd.

Class
    Foam::functionEntries::removeEntry

Description
    Remove a dictionary entry.

    The \c \#remove directive takes a list or a single wordRe.
    For example,
    \verbatim
        #remove entry0
        #remove ( entry1 entry2 entry3 otherEntry )
        #remove "entry[1-3]"
        #remove ( "entry[1-3]" otherEntry )
    \endverbatim

    The removal only occurs in the current context.
    Removing sub-entries or parent entries is not supported.

SourceFiles
    removeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef removeEntry_H
#define removeEntry_H

#include "db/dictionary/functionEntries/functionEntry/functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class removeEntry Declaration
\*---------------------------------------------------------------------------*/

class removeEntry
:
    public functionEntry
{
public:

    //- Runtime type information
    ClassName("remove");


    // Constructors

        //- Disallow default bitwise copy construction
        removeEntry(const removeEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const removeEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
