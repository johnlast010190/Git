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
    initSwakFunctionObject

Description
    Basis for the other function objects here. Evaluates the most common options

SourceFiles
    initSwakFunctionObject.C

Contributors/Copyright:
    2008-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef initSwakFunctionObject_H
#define initSwakFunctionObject_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "primitives/bools/lists/boolList.H"
#include "fields/Fields/Field/Field.H"

#include "include/swak.H"
#include "helpers/switchableNotImplemented.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class initSwakFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class objectRegistry;

class initSwakFunctionObject
:
    public functionObject
{
    // Private data

        dictionary dict_;


public:

    //- Runtime type information
    TypeName("initSwakFunctionObject");


    // Constructors

        //- Construct from components
        initSwakFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        initSwakFunctionObject(const initSwakFunctionObject&) = delete;


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
        //- Interface for old functionObject
        virtual bool execute();
#endif

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        //- Interface for new functionObject
        virtual bool execute();
        virtual bool write();
#endif

        //- Read and set the function object if it's data has changed
        virtual bool read(const dictionary& dict);

#ifdef FOAM_FUNCTIONOBJECT_NEEDS_UPDATEMESH
        //- Make functionObject happy but fail
        //  because this is not an implementation
        virtual void topoChange(const polyTopoChangeMap&)
        {
            switchableNotImplemented
            (
                "initSwakFunctionObject::topoChange",
                dict_
            );
        }

        //- Make functionObject happy but fail
        //  because this is not an implementation
        virtual void movePoints(const polyMesh&)
        {
            switchableNotImplemented
            (
                "initSwakFunctionObject::movePoints",
                dict_
            );
        }
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const initSwakFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
