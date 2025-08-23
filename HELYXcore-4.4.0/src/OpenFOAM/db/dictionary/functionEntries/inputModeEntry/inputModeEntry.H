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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionEntries::inputModeEntry

Description
    Specify the input mode when reading dictionaries, expects
    a single word to follow.

    An example of \c \#inputMode directive:
    \verbatim
        #inputMode merge
    \endverbatim

    The possible input modes:
      - \par merge
        merge sub-dictionaries when possible
      - \par overwrite
        keep last entry and silently remove previous ones
      - \par protect
        keep initial entry and silently ignore subsequent ones
      - \par warn
        keep initial entry and warn about subsequent ones
      - \par error
        issue a FatalError for duplicate entries
      - \par default
        The default treatment - currently identical to \b merge.

      Note that the clear() method resets to the default mode (merge).

SourceFiles
    inputModeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef inputModeEntry_H
#define inputModeEntry_H

#include "db/dictionary/functionEntries/functionEntry/functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class inputModeEntry Declaration
\*---------------------------------------------------------------------------*/

class inputModeEntry
:
    public functionEntry
{
        //- The input mode options
        enum inputMode
        {
            MERGE,
            OVERWRITE,
            PROTECT,
            WARN,
            ERROR
        };

        //- The current input mode
        static inputMode mode_;


    // Private Member Functions

        //- Read the mode as a word and set enum appropriately
        static void setMode(Istream& is);


public:

    //- Runtime type information
    ClassName("inputMode");


    // Constructors

        //- Disallow default bitwise copy construction
        inputModeEntry(const inputModeEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);

        //- Reset the inputMode to %default (ie, %merge)
        static void clear();

        //- Return true if the inputMode is %merge
        static bool merge();

        //- Return true if the inputMode is %overwrite
        static bool overwrite();

        //- Return true if the inputMode is %protect
        static bool protect();

        //- Return true if the inputMode is %error
        static bool error();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inputModeEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
