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

\*---------------------------------------------------------------------------*/

#include "meshToMesh/meshToMeshAdjustTimeStep/meshToMeshAdjustTimeStepFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meshToMeshAdjustTimeStepFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        meshToMeshAdjustTimeStepFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::
meshToMeshAdjustTimeStepFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    meshToMesh_
    (
        refCast<const fvMeshTopoChangers::meshToMesh>(mesh_.topoChanger())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::
~meshToMeshAdjustTimeStepFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


Foam::scalar
Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::timeToNextWrite()
{
    const scalarList& times = meshToMesh_.times();
    const scalar userTimeValue = time_.timeToUserTime(time_.value());

    if (userTimeValue + meshToMesh_.timeDelta() < times.last())
    {
        forAll(times, i)
        {
            if (times[i] > userTimeValue + meshToMesh_.timeDelta())
            {
                return time_.userTimeToTime(times[i] - userTimeValue);
            }
        }
    }

    return VGREAT;
}


bool Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::execute()
{
    return true;
}


bool Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject::write()
{
    return true;
}


// ************************************************************************* //
