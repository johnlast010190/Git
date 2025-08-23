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
    functionObjectListProxy

Description
    Basis for the other function objects here. Evaluates the most common options

SourceFiles
    functionObjectListProxy.C

Contributors/Copyright:
    2011-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef functionObjectListProxy_H
#define functionObjectListProxy_H

#include "db/functionObjects/functionObjectList/functionObjectList.H"
#include "simpleFunctionObject/simpleFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class functionObjectListProxy Declaration
\*---------------------------------------------------------------------------*/

class functionObjectListProxy
:
    public simpleFunctionObject
{
protected:

    // Protected Data

        //- The actual function Objects that this encapsulates
        autoPtr<functionObjectList> functions_;


    // Protected Member Functions

        //- Get the functionObject list (read if not yet done so)
        functionObjectList& functions();

        //- Initialize the function-objects
        virtual void initFunctions();


public:

    //- Runtime type information
    TypeName("functionObjectListProxy");


    // Constructors

        //- Construct from components
        functionObjectListProxy
        (
            const word&,
            const Time&,
            const dictionary&,
            // to switch off for FO that would require a virtual method
            bool allowReadingDuringConstruction=true
        );

        //- Disallow default bitwise copy construction
        functionObjectListProxy(const functionObjectListProxy&) = delete;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        //- Execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);

        //- Write from the simpleFunctionObjects
        virtual void writeSimple();


    // Member Operators

    //- Disallow default bitwise assignment
    void operator=(const functionObjectListProxy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
