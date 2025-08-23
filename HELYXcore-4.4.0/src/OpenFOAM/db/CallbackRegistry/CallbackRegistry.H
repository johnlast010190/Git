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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::CallbackRegistry

Description
    Base class with which callbacks are registered.

    Derive from this class and extend by adding the appropriate callback
    functions that loop and call the callback functions for each entry.

SourceFiles
    CallbackRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef CallbackRegistry_H
#define CallbackRegistry_H

#include "containers/LinkedLists/user/UIDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CallbackType>
class Callback;


/*---------------------------------------------------------------------------*\
                       Class CallbackRegistry Declaration
\*---------------------------------------------------------------------------*/

TemplateName(CallbackRegistry);


/*---------------------------------------------------------------------------*\
                       Class CallbackRegistry Declaration
\*---------------------------------------------------------------------------*/

template<class CallbackType>
class CallbackRegistry
:
    public CallbackRegistryName,
    public UIDLList<CallbackType>
{
public:

    // Constructors

        //- Construct null
        CallbackRegistry();

        //- Disallow default bitwise copy construction
        CallbackRegistry(const CallbackRegistry&) = delete;


    //- Destructor
    virtual ~CallbackRegistry();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CallbackRegistry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/CallbackRegistry/CallbackRegistry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
