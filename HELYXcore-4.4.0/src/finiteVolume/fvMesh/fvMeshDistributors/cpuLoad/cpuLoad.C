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

\*---------------------------------------------------------------------------*/

#include "cpuLoad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cpuLoad, 0);
    optionalCpuLoad optionalCpuLoad::optionalCpuLoad_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuLoad::cpuLoad(const fvMesh& mesh, const word& name)
:
    volScalarField::Internal
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cpuLoad::~cpuLoad()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::optionalCpuLoad& Foam::optionalCpuLoad::New
(
    const fvMesh& mesh,
    const word& name,
    const bool loadBalancing
)
{
    if (loadBalancing)
    {
        if
        (
            mesh.thisDb().objectRegistry::template
            foundObject<cpuLoad>
            (
                name
            )
        )
        {
            return mesh.thisDb().objectRegistry::template
            lookupObjectRef<cpuLoad>
            (
                name
            );
        }
        else
        {
            if (cpuLoad::debug)
            {
                InfoInFunction
                    << "constructing " << name
                    << " for region " << mesh.name() << endl;
            }

            cpuLoad* cpuLoadPtr(new cpuLoad(mesh, name));

            regIOobject::store(cpuLoadPtr);

            return *cpuLoadPtr;
        }
    }
    else
    {
        return optionalCpuLoad::optionalCpuLoad_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuLoad::reset()
{
    cpuTime_.cpuTimeIncrement();
}


void Foam::cpuLoad::cpuTimeIncrement(const label celli)
{
    operator[](celli) += cpuTime_.cpuTimeIncrement();
}


// ************************************************************************* //

