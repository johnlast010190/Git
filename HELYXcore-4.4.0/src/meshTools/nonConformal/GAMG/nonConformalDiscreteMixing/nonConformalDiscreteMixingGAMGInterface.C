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

#include "nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterface.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersection.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(nonConformalLduInterface);

    defineTypeNameAndDebug(nonConformalDiscreteMixingGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        nonConformalDiscreteMixingGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingGAMGInterface::
nonConformalDiscreteMixingGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface(index, coarseInterfaces),
    fineInterface_(refCast<const nonConformalLduInterface>(fineInterface)),
    intersectionsPtr_(nullptr),
    patchIDOffset_(index - fineInterface_.nbrPatch().nbrPatchID())
{
    // Construct face agglomeration from cell agglomeration
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        forAll(localRestrictAddressing, ffi)
        {
            const label curMaster = localRestrictAddressing[ffi];

            auto fnd = masterToCoarseFace.find(curMaster);
            if (fnd == masterToCoarseFace.end())
            {
                // New coarse face
                const label coarseFacei = dynFaceCells.size();
                dynFaceRestrictAddressing.append(coarseFacei);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseFacei);
            }
            else
            {
                // Already have coarse face
                dynFaceRestrictAddressing.append(fnd());
            }
        }

        faceCells_.transfer(dynFaceCells);
        faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
    }

    // On the owner side, construct the discrete mixing intersection engine.
    if (fineInterface_.owner())
    {
        // Construct also the neighbour side agglomeration (as the neighbour
        // interface would do), repeating the exact operation above using
        // neighbourRestrictAddressing instead of localRestrictAddressing).

        labelList nbrFaceRestrictAddressing;
        {
            // From face to coarse face
            DynamicList<label> dynNbrFaceRestrictAddressing
            (
                neighbourRestrictAddressing.size()
            );

            Map<label> masterToCoarseFace(neighbourRestrictAddressing.size());

            forAll(neighbourRestrictAddressing, ffi)
            {
                const label curMaster = neighbourRestrictAddressing[ffi];

                auto fnd = masterToCoarseFace.find(curMaster);
                if (fnd == masterToCoarseFace.end())
                {
                    // New coarse face
                    const label coarseFacei = masterToCoarseFace.size();
                    dynNbrFaceRestrictAddressing.append(coarseFacei);
                    masterToCoarseFace.insert(curMaster, coarseFacei);
                }
                else
                {
                    // Already have coarse face
                    dynNbrFaceRestrictAddressing.append(fnd());
                }
            }

            nbrFaceRestrictAddressing.transfer(dynNbrFaceRestrictAddressing);
        }

        intersectionsPtr_.reset
        (
            new discreteMixingIntersection
            (
                fineInterface_.intersection(),
                faceRestrictAddressing_,
                nbrFaceRestrictAddressing
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingGAMGInterface::
~nonConformalDiscreteMixingGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField>
Foam::nonConformalDiscreteMixingGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    const nonConformalDiscreteMixingGAMGInterface& nbr =
        dynamic_cast<const nonConformalDiscreteMixingGAMGInterface&>
        (
            nbrPatch()
        );
    const labelUList& nbrFaceCells = nbr.faceCells();

    tmp<labelField> tpnf(new labelField(nbrFaceCells.size()));
    labelField& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }

    return tpnf;
}


// ************************************************************************* //
