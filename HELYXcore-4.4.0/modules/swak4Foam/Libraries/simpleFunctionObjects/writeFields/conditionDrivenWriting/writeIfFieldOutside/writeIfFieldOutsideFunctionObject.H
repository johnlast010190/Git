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
    writeIfFieldOutsideFunctionObject

Description
    Checks whether a field exceeds given limits (min and max). Triggers write
    if that is the case

SourceFiles
    writeIfFieldOutsideFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id: writeIfFieldOutsideFunctionObject.H,v 78b0d113b99b 2013-02-25 16:12:41Z bgschaid $
\*---------------------------------------------------------------------------*/

#ifndef writeIfFieldOutsideFunctionObject_H
#define writeIfFieldOutsideFunctionObject_H

#include "writeFields/conditionDrivenWriting/conditionDrivenWritingFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class writeIfFieldOutsideFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeIfFieldOutsideFunctionObject
:
    public conditionDrivenWritingFunctionObject
{
    // Private Data

        //- Name of the field that is to be checked
        word fieldName_;

        //- Maximum value that is to be looked for
        scalar maximum_;

        //- Minimum value that is to be looked for
        scalar minimum_;


    // Private Member Functions

        //- Do the actual checking
        template <class T>
        scalar getMin() const;

        //- Do the actual checking
        template <class T>
        scalar getMax() const;

        //- Wrapper for the checking
        template <class T>
        bool check() const;


protected:

    // Protected Member Functions

        //- Should writing start now?
        virtual bool checkStartWriting();


public:

    //- Runtime type information
    TypeName("writeIfFieldOutside");


    // Constructors

        //- Construct from components
        writeIfFieldOutsideFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeIfFieldOutsideFunctionObject
        (
            const writeIfFieldOutsideFunctionObject&
        ) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeIfFieldOutsideFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "writeIfFieldOutsideFunctionObjectTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
