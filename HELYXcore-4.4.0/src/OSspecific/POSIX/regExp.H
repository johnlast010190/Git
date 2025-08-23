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
    Foam::regExp

Description
    Wrapper around POSIX extended regular expressions.

    The PCRE '(?i)' extension is provided to compile the regular expression
    as being case-insensitive.

See also
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from \c Perl and \c sed regular expressions.

SourceFiles
    regExp.C

\*---------------------------------------------------------------------------*/

#ifndef regExp_H
#define regExp_H

#include <regex.h>
#include <string>
#include "containers/Lists/List/ListDecl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regExp Declaration
\*---------------------------------------------------------------------------*/

class regExp
{
    // Private data

        //- Precompiled regular expression
        regex_t* preg_;


    // Private Member Functions

        //- Return true if it matches and sets the sub-groups matched.
        bool matchGrouping
        (
            const std::string& text,
            List<std::string>& groups
        ) const;


public:

    // Static Member Functions

        //- Is character a regular expression meta-character?
        //  any character: '.' \n
        //  quantifiers: '*', '+', '?' \n
        //  grouping: '(', '|', ')' \n
        //  range: '[', ']' \n
        //
        //  Don't bother checking for '{digit}' bounds
        inline static bool meta(const char c);


    // Constructors

        //- Construct null
        inline regExp();

        //- Construct from character array
        inline explicit regExp(const char* pattern);

        //- Construct from string
        inline explicit regExp(const std::string& pattern);

        //- Construct from character array, optionally ignore case
        inline regExp(const char* pattern, bool ignoreCase);

        //- Construct from string, optionally ignore case
        inline regExp(const std::string& pattern, bool ignoreCase);

        //- Disallow default bitwise copy construction
        regExp(const regExp&) = delete;


    //- Destructor
    inline ~regExp();


    // Member functions

      // Access

        //- Return true if a precompiled expression does not exist
        inline bool empty() const;

        //- Does a precompiled expression exist?
        inline bool exists() const;

        //- The number of capture groups for a non-empty expression
        inline unsigned ngroups() const;


      // Editing

        //- Compile pattern into a regular expression, optionally ignore case.
        //  \return True if the pattern was compiled
        bool set(const char* pattern, bool ignoreCase=false);

        //- Compile pattern into a regular expression, optionally ignore case
        //  \return True if the pattern was compiled
        bool set(const std::string& pattern, bool ignoreCase=false);

        //- Clear expression.
        //  \return True if expression had existed prior to the clear.
        bool clear();


      // Matching/Searching

        //- Find position within the text.
        //  \return The index where it begins or string::npos if not found
        std::string::size_type find(const std::string& text) const;

        //- True if the regex matches the entire text.
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        bool match(const std::string& text) const;

        //- True if the regex matches the text, set the sub-groups matched.
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        bool match(const std::string& text, List<std::string>& groups) const;

        //- Return true if the regex was found within the text
        inline bool search(const std::string& text) const;


    // Member Operators

        //- Perform match on text
        inline bool operator()(const std::string& text) const;

        //- Assign and compile pattern from a character array
        //  Always case sensitive
        inline void operator=(const char* pattern);

        //- Assign and compile pattern from string
        //  Always case sensitive
        inline void operator=(const std::string& pattern);

        //- Disallow default bitwise assignment
        void operator=(const regExp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regExpI.H"

#endif

// ************************************************************************* //
