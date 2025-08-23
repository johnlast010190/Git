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
    executeIfWriteTimeFunctionObject

Description
    This function objects only executes its children only if the current time
    is a writeTime

SourceFiles
    executeIfWriteTimeFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfWriteTimeFunctionObject_H
#define executeIfWriteTimeFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"
#include "regExp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class executeIfWriteTimeFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfWriteTimeFunctionObject
:
    public conditionalFunctionObjectListProxy
{
protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfWriteTime");


    // Constructors

        //- Construct from components
        executeIfWriteTimeFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfWriteTimeFunctionObject
        (
            const executeIfWriteTimeFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfWriteTimeFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
