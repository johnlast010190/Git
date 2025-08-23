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

Class
    Foam::ensight::VarName

Description
    Specification of a valid Ensight variable-name.

    \code
        space !#%()*+,-./;@[]^
    \endcode
    If spaces exist, they must be quoted

    Variable names cannot start with a digit.
    The maximum part name length (in GUI: 49 chars) - not addresssed.
    The maximum variable length (in GUI: 49 chars) - not addresssed.

\*---------------------------------------------------------------------------*/

#ifndef ensightVarName_H
#define ensightVarName_H

#include "primitives/strings/word/word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ensight
{

/*---------------------------------------------------------------------------*\
                      Class ensight::VarName Declaration
\*---------------------------------------------------------------------------*/

class VarName
:
    public word
{
    // Private Member Functions

        //- Strip invalid characters
        inline void stripInvalid();

public:

    // Constructors

        //- Construct as copy
        inline explicit VarName(const VarName&);

        //- Construct as copy of character array
        inline explicit VarName(const char*);

        //- Construct as copy of std::string
        inline explicit VarName(const std::string&);


    // Member functions

        //- Is this character valid for an ensight var-name
        inline static bool valid(char);


    // Member operators

    // Assignment (disabled)
        void operator=(const word&) = delete;
        void operator=(const string&) = delete;
        void operator=(const std::string&) = delete;
        void operator=(const char*) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensight
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensight/name/ensightVarNameI.H"

#endif

// ************************************************************************* //
