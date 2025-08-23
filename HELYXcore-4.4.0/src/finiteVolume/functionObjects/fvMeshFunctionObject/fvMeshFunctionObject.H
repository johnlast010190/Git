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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::fvMeshFunctionObject

Description
    Specialization of Foam::functionObject for an Foam::fvMesh, providing a
    reference to the Foam::fvMesh.

    If the selected region is not an Foam::fvMesh a Foam::FatalError will be
    generated.

See also
    Foam::regionFunctionObject
    Foam::functionObject

SourceFiles
    fvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fvMeshFunctionObject_H
#define functionObjects_fvMeshFunctionObject_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fvMeshFunctionObject
:
    public regionFunctionObject
{
protected:

    // Protected member data

        //- Reference to the fvMesh
        const fvMesh& mesh_;


private:

    // Private Member Functions

        //- Select the mesh to store in member data
        const fvMesh& selectMesh
        (
            const Time&,
            const objectRegistry&
        ) const;


public:

    //- Runtime type information
    TypeName("fvMeshFunctionObject");

    //- Return the base directory for output
    fileName outputFileDir() const;

    // Constructors

        //- Construct from Time and dictionary
        fvMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        fvMeshFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fvMeshFunctionObject(const fvMeshFunctionObject&) = delete;


    //- Destructor
    virtual ~fvMeshFunctionObject();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
