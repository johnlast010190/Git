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
    (c) 2019-2024 Engys Ltd.

Authors
    2019. Nikolaos Magoulas (Engys Ltd.). All rights reserved.

Class
    Foam::dynamicSwitchingFvMesh

Description
    The dynamicSwitchingFvMesh allows for a mesh to define multiple boundary
    types per boundary patch and switch between each type depending on the
    current state.

SourceFiles
    dynamicSwitchingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicSwitchingFvMesh_H
#define dynamicSwitchingFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dynamicSwitchingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicSwitchingFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Name of the currently activated database
        word currentRegistryName_;

        //- List of the registries' names
        wordList registries_;

        //- A list of databases that holds a list of the state of each patch
        //  Registries are sorted in the same way as in the registries wordList
        boolListList regPatchState_;


    // Private Member Functions

        //- Read the registries
        wordList readRegistries(const IOobject& io);

        //- Read the states of each patch
        boolListList readPatches
        (
            const IOobject& io,
            const wordList& registries
        );


public:

    //- Runtime type information
    TypeName("dynamicSwitchingFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicSwitchingFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicSwitchingFvMesh(const dynamicSwitchingFvMesh&) = delete;


    //- Destructor
    ~dynamicSwitchingFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update the mesh to the provided object registry
        virtual bool update(const word& registryName);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicSwitchingFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
