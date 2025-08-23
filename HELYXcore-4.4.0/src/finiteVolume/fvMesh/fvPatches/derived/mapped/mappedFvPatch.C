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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, mappedFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::mappedFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return this->patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::mappedFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    // This is the region-segregated case - region-coupled is below in
    // regionCoupledInternalFieldTransfer
    if (isInterface())
    {
        if (sampleRegion() == this->boundaryMesh().mesh().name())
        {
            // Patch is coupled back to another patch in same mesh
            return nbrPatch().patchInternalField(iF);
        }
        else
        {
            // Patch is coupled but the field being solved is not.
            // Return a dummy array
            return tmp<labelField>(new labelField(nbrPatch().size(), -1));
        }
    }
    else
    {
        // Coupling BC should not be used for this type of sampling, or
        // region-coupled boundary is being called by segregated solver where
        // neighbour mesh does not exist.
        // Return dummy array
        return tmp<labelField>(new labelField(this->size(), -1));
    }

}


Foam::tmp<Foam::labelField>
Foam::mappedFvPatch::regionCoupledInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    HELYX_ASSERT(isInterface())
    {
        FatalErrorInFunction
            << "Region-coupled interface expected"
            << exit(FatalError);
    }

    return nbrPatch().patchInternalField(iF);
}

// ************************************************************************* //
