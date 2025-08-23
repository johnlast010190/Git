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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "motionInterpolation/patchCorrected/patchCorrectedInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchCorrectedInterpolation, 0);

    addToRunTimeSelectionTable
    (
        motionInterpolation,
        patchCorrectedInterpolation,
        Istream
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::patchCorrectedInterpolation::getPatchGroups
(
    Istream& entry
) const
{
    List<List<word>> patchGroupNames(entry);

    labelListList patchGroups(patchGroupNames.size());

    forAll(patchGroupNames, patchI)
    {
        patchGroups[patchI].resize(patchGroupNames[patchI].size());

        forAll(patchGroupNames[patchI], patchJ)
        {
            patchGroups[patchI][patchJ] =
                mesh().boundaryMesh().findPatchID
                (
                    patchGroupNames[patchI][patchJ]
                );

            if (patchGroups[patchI][patchJ] == -1)
            {
                FatalErrorInFunction
                    << "patch \"" << patchGroupNames[patchI][patchJ]
                    << "\" not found" << exit(FatalError);
            }
        }
    }

    return patchGroups;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchCorrectedInterpolation::patchCorrectedInterpolation
(
    const fvMesh& mesh,
    Istream& entry
)
:
    motionInterpolation(mesh, entry),
    patchGroups_(getPatchGroups(entry))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchCorrectedInterpolation::~patchCorrectedInterpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::patchCorrectedInterpolation::interpolate
(
    const volScalarField& cellDisplacement,
    pointScalarField& pointDisplacement
) const
{
    interpolateType(cellDisplacement, pointDisplacement);
}


void Foam::patchCorrectedInterpolation::interpolate
(
    const volVectorField& cellDisplacement,
    pointVectorField& pointDisplacement
) const
{
    interpolateType(cellDisplacement, pointDisplacement);
}


// ************************************************************************* //
