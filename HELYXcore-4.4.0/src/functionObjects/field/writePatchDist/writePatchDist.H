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
    (c) 2021-2024 Engys Ltd.
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::writePatchDist

Group
    grpForcesFunctionObjects

Description
    Writes the patch distance field

    All wall patch patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    writePatchDist1
    {
        type        writePatchDist;
        patches     (ground);
        executeControl writeTime;
        executeInterval 1;
        writeControl outputTimeAndEnd;
        writeInterval 1;
        functionObjectLibs ( "libfieldFunctionObjects.so" );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: writePatchDist    | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    writePatchDist.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writePatchDist_H
#define functionObjects_writePatchDist_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fvMesh/fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writePatchDist Declaration
\*---------------------------------------------------------------------------*/

class writePatchDist
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        const fvMesh& mesh_;

        //- Optional list of patches to process
        labelHashSet patchSet_;


public:

    //- Runtime type information
    TypeName("writePatchDist");


    // Constructors

        //- Construct from Time and dictionary
        writePatchDist
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writePatchDist(const writePatchDist&) = delete;


    //- Destructor
    virtual ~writePatchDist();


    // Member Functions

        //- Read the writePatchDist data
        virtual bool read(const dictionary&);

        //- Calculate the patch distance
        virtual bool execute();

        //- Write the patch distance
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writePatchDist&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
