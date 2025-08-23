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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/Time/Time.H"
#include "ensight/output/ensightOutput.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
int Foam::functionObjects::ensightWrite::writeVolField
(
    const fvMeshSubset& proxy,
    const word& inputName,
    int& state
)
{
    // State: return 0 (not-processed), -1 (skip), +1 ok
    typedef VolField<Type> VolFieldType;
    const fvMesh& baseMesh = proxy.baseMesh();

    // Already done, or not available
    if (state || !baseMesh.foundObject<VolFieldType>(inputName))
    {
        return state;
    }

    const auto* fieldptr = &baseMesh.lookupObject<VolFieldType>(inputName);
    auto tfield = fvMeshSubsetProxy::interpolate(proxy, *fieldptr);
    const auto& field = tfield();

    autoPtr<ensightFile> os = ensCase().newData<Type>(inputName);
    ensightOutput::writeField<Type>
    (
        field,
        ensMesh(),
        os
    );

    Log << " " << inputName;

    state = +1;
    return state;
}


// ************************************************************************* //
