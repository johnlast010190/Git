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
    executeIfStartTimeFunctionObject

Description
    This function objects only executes its children only if the current time
    is the startTime

SourceFiles
    executeIfStartTimeFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfStartTimeFunctionObject_H
#define executeIfStartTimeFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"
#include "regExp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class executeIfStartTimeFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfStartTimeFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- Execute only if this is the start time
        bool runIfStartTime_;

        //- If this is a restart execute as well
        bool executeOnRestart_;

        //- The time index at which this was read
        label timeIndexRead_;


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfStartTime");


    // Constructors

        //- Construct from components
        executeIfStartTimeFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfStartTimeFunctionObject
        (
            const executeIfStartTimeFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfStartTimeFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
