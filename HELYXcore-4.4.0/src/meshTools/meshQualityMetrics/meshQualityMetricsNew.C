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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshQualityMetrics/meshQualityMetrics.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "db/IOobjects/IOdictionary/localIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::meshQualityMetrics& Foam::meshQualityMetrics::New
(
    const polyMesh& mesh,
    const word& mqName
)
{
    if (!mesh.foundObject<meshQualityMetrics>(mqName))
    {
        const word dictName = "meshObjects";

        autoPtr<IOobject> meshObjectIO
        (
            new IOobject
            (
                dictName,
                mesh.time().caseSystem(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if (!meshObjectIO->typeHeaderOk<localIOdictionary>())
        {
            // Fallback: look for definition at the top level
            autoPtr<IOobject> defaultRegionIO
            (
                new IOobject
                (
                    dictName,
                    mesh.time().caseSystem(),
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            if (defaultRegionIO->typeHeaderOk<localIOdictionary>())
            {
                meshObjectIO = defaultRegionIO;
            }
        }

        // Use objectPath rather than filePath so a useful error is given if
        // file is not found
        IFstream is(meshObjectIO->objectPath());
        dictionary meshObjDict(is);
        dictionary mqDict = meshObjDict.subDict(mqName);

        autoPtr<meshQualityMetrics> mqPtr
        (
            new meshQualityMetrics(mesh, mqDict, mqName)
        );

        mqPtr.ptr()->store();
    }

    return mesh.lookupObjectRef<meshQualityMetrics>(mqName);
}



// ************************************************************************* //
