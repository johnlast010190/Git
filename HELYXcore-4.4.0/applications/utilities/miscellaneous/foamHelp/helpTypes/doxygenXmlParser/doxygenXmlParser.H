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
    (c) 2012-2014 OpenFOAM Foundation

Class
    Foam::doxygenXmlParser

Description
    Parser for doxygen XML

SourceFiles
    doxygenXmlParser.C

\*---------------------------------------------------------------------------*/

#ifndef doxygenXmlParser_H
#define doxygenXmlParser_H

#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class doxygenXmlParser Declaration
\*---------------------------------------------------------------------------*/

class doxygenXmlParser
:
    public dictionary
{

public:

    //- Construct from components
    doxygenXmlParser
    (
        const fileName& fName,
        const string& startTag,
        const string& searchStr,
        const bool exactMatch,
        const word& ext
    );


    // Member functions

        //- Skip past a block
        void skipBlock(IFstream& is, const word& blockName) const;

        //- Skip forward to block
        void skipForward(IFstream& is, const word& blockName) const;

        //- Return the entry
        template<class Type>
        void getEntry(IFstream& is, Type& entry) const;

        //- Return the entry value
        template<class Type>
        void getValue(IFstream& is, Type& entry) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "doxygenXmlParserTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
