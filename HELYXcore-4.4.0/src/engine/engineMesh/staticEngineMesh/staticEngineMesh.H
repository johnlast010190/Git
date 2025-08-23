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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::staticEngineMesh

Description
    Foam::staticEngineMesh

SourceFiles
    staticEngineMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticEngineMesh_H
#define staticEngineMesh_H

#include "engineMesh/engineMesh/engineMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class staticEngineMesh Declaration
\*---------------------------------------------------------------------------*/

class staticEngineMesh
:
    public engineMesh
{
public:

    //- Runtime type information
    TypeName("static");


    // Constructors

        //- Construct from IOobject
        staticEngineMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        staticEngineMesh(const staticEngineMesh&) = delete;


    //- Destructor
    ~staticEngineMesh();


    // Member Functions

        //- Move mesh
        virtual void moveMesh() override
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const staticEngineMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
