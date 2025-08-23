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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::functionEntries::includeFuncEntry

Description
    Specify a functionObject dictionary file to include, expects the
    functionObject name to follow with option arguments (without quotes).

    Searches for functionObject dictionary file in user/group/shipped
    directories allowing for version-specific and version-independent files
    using the following hierarchy:
    - \b user settings:
      - ~/.HELYXcore/\<VERSION\>/caseDicts/postProcessing
      - ~/.HELYXcore/caseDicts/postProcessing
    - \b group (site) settings (when $HELYX_CUSTOM_DATA_DIR is set):
      - $HELYX_CUSTOM_DATA_DIR/\<VERSION\>/caseDicts/postProcessing
      - $HELYX_CUSTOM_DATA_DIR/caseDicts/postProcessing
    - \b group (site) settings (when $HELYX_CUSTOM_DATA_DIR is not set):
      - $HELYX_PROJECT_DIR/../site/\<VERSION\>/caseDicts/postProcessing
      - $HELYX_PROJECT_DIR/../site/caseDicts/postProcessing
    - \b other (project) settings:
      - $HELYX_PROJECT_DIR/etc/caseDicts/postProcessing

    The optional field arguments included in the name are inserted in 'field' or
    'fields' entries in the functionObject dictionary and included in the name
    of the functionObject entry to avoid conflict.

    Examples:
    \verbatim
        #includeFunc Q
        #includeFunc components(U)
        #includeFunc mag(Ux)
        #includeFunc mag(p)
    \endverbatim

See also
    Foam::functionObjectList

SourceFiles
    includeFuncEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeFuncEntry_H
#define includeFuncEntry_H

#include "db/dictionary/functionEntries/functionEntry/functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeFuncEntry Declaration
\*---------------------------------------------------------------------------*/

class includeFuncEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    ClassName("includeFunc");


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
