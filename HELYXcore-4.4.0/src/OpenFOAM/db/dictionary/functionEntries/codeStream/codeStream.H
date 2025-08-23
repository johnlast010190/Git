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
    (c) 2011-2020 OpenFOAM Foundation

Class
    Foam::functionEntries::codeStream

Description
    Dictionary entry that contains C++ OpenFOAM code that is compiled to
    generate the entry itself. So
    - codeStream reads three entries: 'code', 'codeInclude' (optional),
    'codeOptions' (optional)
    and uses those to generate library sources inside \c codeStream/
    - these get compiled using 'wmake libso'
    - the resulting library is loaded in executed with as arguments
    \code
        (const dictionary& dict, Ostream& os)
    \endcode
    where the dictionary is the current dictionary.
    - the code has to write into Ostream which is then used to construct
    the actual dictionary entry.

    E.g. to set the internal field of a field:

    \verbatim
    internalField  #codeStream
    {
        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            scalarField fld(mesh.nCells(), 12.34);
            fld.writeEntry("", os);
        #};

        //- Optional:
        codeInclude
        #{
            #include "cfdTools/general/include/fvCFD.H"
        #};

        //- Optional:
        codeOptions
        #{
            -I$(HELYX_PROJECT_DIR)/src/finiteVolume
        #};
    };
    \endverbatim

    Note the \c \#{ ... \c \#} syntax is a 'verbatim' input mode that allows
    inputting strings with embedded newlines.

    Limitations:
    - '~' symbol not allowed inside the code sections.
    - probably some other limitations (uses string::expand which expands
      \c \$ and \c ~ sequences)

Note
    The code to be compiled is stored under the local \c codeStream directory
    with a subdirectory name corresponding to the SHA1 of the contents.

    The corresponding library code is located under the local
    \c codeStream/platforms/$HELYX_OPTIONS/lib directory in a library
    \c libcodeStream_SHA1.so

SourceFiles
    codeStream.C

\*---------------------------------------------------------------------------*/

#ifndef codeStream_H
#define codeStream_H

#include "db/dictionary/functionEntries/functionEntry/functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend classes
class dlLibraryTable;

namespace functionEntries
{

class calcEntry;

/*---------------------------------------------------------------------------*\
                          Class codeStream Declaration
\*---------------------------------------------------------------------------*/

class codeStream
:
    public functionEntry
{

    //- Interpreter function type
    typedef void (*streamingFunctionType)(Ostream&, const dictionary&);


    // Private Member Functions

        //- Return true for master only reading of global dictionary
        static bool masterOnlyRead(const dictionary& dict);

        //- Construct, compile, load and return streaming function
        static streamingFunctionType getFunction
        (
            const dictionary& parentDict,
            const dictionary& codeDict
        );


public:

    // Static data members

        //- Name of the C code template to be used
        static const word codeTemplateC;


    // Related types

        //- Declare friendship with the calcEntry class
        friend class calcEntry;


    //- Runtime type information
    ClassName("codeStream");


    // Constructors

        //- Disallow default bitwise copy construction
        codeStream(const codeStream&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const dictionary& parentDict,
            primitiveEntry& entry,
            Istream& is
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codeStream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
