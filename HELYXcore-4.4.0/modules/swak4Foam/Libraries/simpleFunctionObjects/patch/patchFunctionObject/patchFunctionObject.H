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
    patchFunctionObject

Description
    Basis for all the function objects working on "per patch" basis

SourceFiles
    patchFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef patchFunctionObject_H
#define patchFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class patchFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchFunctionObject
:
    public timelineFunctionObject
{
protected:

    // Protected Data

        //- Names of the patches to use
        wordList patchNames_;

        //- indizes of the patches
        List<label> patchIndizes_;


    // Protected Member Functions

        virtual stringList columnNames();


public:

    //- Runtime type information
    TypeName("patchFunctionObject");


    // Constructors

        //- Construct from components
        patchFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchFunctionObject(const patchFunctionObject&) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
