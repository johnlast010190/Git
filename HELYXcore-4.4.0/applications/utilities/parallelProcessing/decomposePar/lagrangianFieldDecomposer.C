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
    (c) 2011-2016 OpenFOAM Foundation

Description
    Lagrangian field decomposer.

\*---------------------------------------------------------------------------*/

#include "lagrangianFieldDecomposer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianFieldDecomposer::lagrangianFieldDecomposer
(
    const polyMesh& mesh,
    const polyMesh& procMesh,
    const labelList& faceProcAddressing,
    const labelList& cellProcAddressing,
    const word& cloudName,
    const Cloud<indexedParticle>& lagrangianPositions,
    const List<SLList<indexedParticle*>*>& cellParticles
)
:
    procMesh_(procMesh),
    positions_(procMesh, cloudName, IDLList<passiveParticle>()),
    particleIndices_(lagrangianPositions.size())
{
    label pi = 0;

    // faceProcAddressing not required currently
    // labelList decodedProcFaceAddressing(faceProcAddressing.size());

    // forAll(faceProcAddressing, i)
    // {
    //     decodedProcFaceAddressing[i] = mag(faceProcAddressing[i]) - 1;
    // }

    forAll(cellProcAddressing, procCelli)
    {
        label celli = cellProcAddressing[procCelli];

        if (cellParticles[celli])
        {
            SLList<indexedParticle*>& particlePtrs = *cellParticles[celli];

            forAllConstIter(SLList<indexedParticle*>, particlePtrs, iter)
            {
                const indexedParticle& ppi = *iter();
                particleIndices_[pi++] = ppi.index();

                // label mappedTetFace = findIndex
                // (
                //     decodedProcFaceAddressing,
                //     ppi.tetFace()
                // );

                // if (mappedTetFace == -1)
                // {
                //     FatalErrorInFunction
                //         << "Face lookup failure." << nl
                //         << abort(FatalError);
                // }

                positions_.append
                (
                    new passiveParticle
                    (
                        procMesh,
                        ppi.position(),
                        procCelli,
                        false
                    )
                );
            }
        }
    }

    particleIndices_.setSize(pi);

    IOPosition<Cloud<passiveParticle>>(positions_).write();
}


// ************************************************************************* //
