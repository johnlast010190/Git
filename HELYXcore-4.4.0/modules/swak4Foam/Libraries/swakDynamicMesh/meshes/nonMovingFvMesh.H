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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::nonMovingFvMesh

Description
    Foam::nonMovingFvMesh - Basically the staticFvMesh-class but it pretends
    to have moved

SourceFiles
    nonMovingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef nonMovingFvMesh_H
#define nonMovingFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class nonMovingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class nonMovingFvMesh
:
    public dynamicFvMesh
{
public:

    //- Runtime type information
    TypeName("nonMovingFvMesh");


    // Constructors

        //- Construct from IOobject
        nonMovingFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        nonMovingFvMesh(const nonMovingFvMesh&) = delete;


    //- Destructor
    ~nonMovingFvMesh();


    // Member Functions

        //- Dummy update function which does not change the mesh
        virtual bool update();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonMovingFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
