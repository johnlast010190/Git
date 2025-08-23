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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd.

Class
    printDictionary

Group
    grpReportFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call printDictionary.write() to calculate&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    printDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_printDictionary_H
#define functionObjects_printDictionary_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class printDictionary Declaration
\*---------------------------------------------------------------------------*/

class printDictionary
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("printDictionary");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        printDictionary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        printDictionary(const printDictionary&) = delete;


    //- Destructor
    virtual ~printDictionary();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const printDictionary&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
