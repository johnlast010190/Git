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
    executeIfParallelSerialFunctionObject

Description
    This function objects only executes its children only if the run is
    serial/parallel

SourceFiles
    executeIfParallelSerialFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfParallelSerialFunctionObject_H
#define executeIfParallelSerialFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"
#include "regExp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class executeIfParallelSerialFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfParallelSerialFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- Execute only if this is a parallel run
        bool runIfParallel_;


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfParallelSerial");


    // Constructors

        //- Construct from components
        executeIfParallelSerialFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfParallelSerialFunctionObject
        (
            const executeIfParallelSerialFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfParallelSerialFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
