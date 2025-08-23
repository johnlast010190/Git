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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::lduRegisteredPrimitiveMesh

Description
    lduPrimitiveMesh which stores an object registry (needed by GAMG).

SourceFiles
    lduRegisteredPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduRegisteredPrimitiveMesh_H
#define lduRegisteredPrimitiveMesh_H

#include "meshes/lduMesh/lduPrimitiveMesh.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduRegisteredPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduRegisteredPrimitiveMesh
:
    public lduPrimitiveMesh
{
    // Private data

        //- Object registry returned by thisDb
        const objectRegistry& obr_;


public:

    // Static data

        //- Declare name of the class and its debug switch
        ClassName("lduRegisteredPrimitiveMesh");


    // Constructors

        //- Construct from components and re-use storage.
        lduRegisteredPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            PtrList<const lduInterface>& primitiveInterfaces,
            const lduSchedule& ps,
            const label comm,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        lduRegisteredPrimitiveMesh(const lduRegisteredPrimitiveMesh&) = delete;


    //- Destructor
    virtual ~lduRegisteredPrimitiveMesh()
    {}


    // Member Functions

        //- Return the object registry
        virtual const objectRegistry& thisDb() const
        {
            return obr_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lduRegisteredPrimitiveMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
