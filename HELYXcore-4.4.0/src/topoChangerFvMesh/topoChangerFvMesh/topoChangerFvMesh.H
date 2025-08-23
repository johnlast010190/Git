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
    Foam::topoChangerFvMesh

Description
    Abstract base class for a topology changing fvMesh

SourceFiles
    topoChangerFvMesh.C
    newTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef topoChangerFvMesh_H
#define topoChangerFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class topoChangerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class topoChangerFvMesh
:
    public dynamicFvMesh
{
protected:

    // Protected Data

        polyTopoChanger topoChanger_;


public:

    //- Runtime type information
    TypeName("topoChangerFvMesh");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit topoChangerFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        topoChangerFvMesh(const topoChangerFvMesh&) = delete;


    //- Destructor
    virtual ~topoChangerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const topoChangerFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
