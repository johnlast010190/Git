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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::removeRegisteredObject

Group
    grpUtilitiesFunctionObjects

Description
    Removes registered objects if present in the database.

Usage
    Example of function object specification:
    \verbatim
    removeRegisteredObject1
    {
        type        removeRegisteredObject;
        libs        ("libutilityFunctionObjects.so");
        ...
        objects     (obj1 obj2);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: removeRegisteredObject | yes |
        objects      | objects to remove       | yes         |
    \endtable

See also
    Foam::functionObject

SourceFiles
    removeRegisteredObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_removeRegisteredObject_H
#define functionObjects_removeRegisteredObject_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "primitives/strings/lists/wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class removeRegisteredObject Declaration
\*---------------------------------------------------------------------------*/

class removeRegisteredObject
:
    public regionFunctionObject
{
    // Private data

        //- Names of objects to control
        wordList objectNames_;


public:

    //- Runtime type information
    TypeName("removeRegisteredObject");


    // Constructors

        //- Construct from Time and dictionary
        removeRegisteredObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        removeRegisteredObject(const removeRegisteredObject&) = delete;


    //- Destructor
    virtual ~removeRegisteredObject();


    // Member Functions

        //- Read the removeRegisteredObject data
        virtual bool read(const dictionary&);

        //- Remove the registered objects
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const removeRegisteredObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
