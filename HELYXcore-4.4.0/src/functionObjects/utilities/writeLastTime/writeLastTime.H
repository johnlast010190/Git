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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::writeLastTime

Description

    Write the database at the final iteration (if needed)

    Example of function object specification:

    \verbatim
    writeLastTime1
    {
        type        writeLastTime;
        functionObjectLibs ("libutilityFunctionObjects.so");
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: writeLastTime | yes |
    \endtable

SeeAlso
    Foam::functionObject

SourceFiles
    writeLastTime.C

\*---------------------------------------------------------------------------*/

#ifndef writeLastTime_H
#define writeLastTime_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "writeLastTime/writeLastTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class writeLastTime Declaration
\*---------------------------------------------------------------------------*/

class writeLastTime
:
    public functionObject
{
    // Private data

        //- Refererence to Db
        const objectRegistry& obr_;


public:

    //- Runtime type information
    TypeName("writeLastTime");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        writeLastTime
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeLastTime(const writeLastTime&) = delete;


    //- Destructor
    virtual ~writeLastTime();


    // Member Functions

        //- Read the writeLastTime data
        virtual bool read(const dictionary&)
        {
            return false;
        }

        //- Execute, currently does nothing
        virtual bool execute()
        {
            return false;
        }

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the writeLastTime
        virtual bool write()
        {
            return false;
        }

        //- Update for changes of mesh
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeLastTime&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
