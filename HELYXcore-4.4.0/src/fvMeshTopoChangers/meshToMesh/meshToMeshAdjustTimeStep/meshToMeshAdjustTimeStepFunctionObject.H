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
    (c) 2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject

Description
    Adjusts time-step for meshToMesh mapping.

SourceFiles
    meshToMeshAdjustTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_meshToMeshAdjustTimeStepFunctionObject_H
#define functionObjects_meshToMeshAdjustTimeStepFunctionObject_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "meshToMesh/fvMeshTopoChangersMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class meshToMeshAdjustTimeStepFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class meshToMeshAdjustTimeStepFunctionObject
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Reference to the meshToMesh fvMeshTopoChanger
        const fvMeshTopoChangers::meshToMesh& meshToMesh_;


public:

    //- Runtime type information
    TypeName("meshToMeshAdjustTimeStep");


    // Constructors

        //- Construct from components
        meshToMeshAdjustTimeStepFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        meshToMeshAdjustTimeStepFunctionObject
        (
            const meshToMeshAdjustTimeStepFunctionObject&
        ) = delete;


    // Destructor
    virtual ~meshToMeshAdjustTimeStepFunctionObject();


    // Member Functions

        //-  Read and reset the timeStep Function1
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Return the time to the next write
        virtual scalar timeToNextWrite();

        //- Do nothing
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMeshAdjustTimeStepFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
