/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2024 Engys Ltd.

Class
    writeFieldsGeneralFunctionObject

Description
    General Class for writting Fields

SourceFiles
    writeFieldsGeneralFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2014 David Huckaby <e.david.huckaby@netl.doe.gov>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeFieldsGeneralFunctionObject_H
#define writeFieldsGeneralFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class writeFieldsGeneralFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeFieldsGeneralFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        //- Name of the fields that are to be written
        wordList fieldNames_;

        //- Name of the fields that are to be written
        wordList cloudNames_;


    // Private Member Functions

        //- Wrapper for the writing
        template <class T>
        bool writeField(const word &name) const;

        //- Wrapper for writing clouds
        template <class T>
        bool writeCloud(const word &name) const;


public:

    // Constructors

        //- Construct from components
        writeFieldsGeneralFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeFieldsGeneralFunctionObject
        (
            const writeFieldsGeneralFunctionObject&
        ) = delete;


    // Member Functions

        bool start();

        void writeSimple();

        const wordList& fieldNames() const
        {
            return fieldNames_;
        }

        const wordList& cloudNames() const
        {
            return cloudNames_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeFieldsGeneralFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "writeFieldsGeneralFunctionObjectTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
