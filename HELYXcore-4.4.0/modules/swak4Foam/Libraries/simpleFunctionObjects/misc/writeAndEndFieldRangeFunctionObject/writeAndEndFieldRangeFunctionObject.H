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
    writeAndEndFieldRangeFunctionObject

Description
    Checks whether a field exceeds given limits (min and max). If the
    limits are exceeded everything is written to disk and the run ends gracefully

SourceFiles
    writeAndEndFieldRangeFunctionObject.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeAndEndFieldRangeFunctionObject_H
#define writeAndEndFieldRangeFunctionObject_H

#include "misc/writeAndEndFunctionObject/writeAndEndFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class writeAndEndFieldRangeFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeAndEndFieldRangeFunctionObject
:
    public writeAndEndFunctionObject
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

        //- Do the actual checking
        virtual bool endRunNow();


public:

    //- Runtime type information
    TypeName("writeAndEndFieldRange");


    // Constructors

        //- Construct from components
        writeAndEndFieldRangeFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeAndEndFieldRangeFunctionObject
        (
            const writeAndEndFieldRangeFunctionObject&
        ) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeAndEndFieldRangeFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "writeAndEndFieldRangeFunctionObjectI.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
