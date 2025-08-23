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

Class
    Foam::limitVolumes

Description

SourceFiles
    limitVolumes.H

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/limitVolumes/limitVolumes.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::limitVolumes
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    if (dict.found("limitVolumeCoeffs"))
    {
        const dictionary& volDict = dict.subDict("limitVolumeCoeffs");

        scalar threshold =
            volDict.lookupOrDefault<scalar>("neiVolFractionThreshold", 2e-2);
        Switch printInformation =
            volDict.lookupOrDefault<Switch>("printInformation", false);

        tmp<volScalarField> tVols
        (
            new volScalarField
            (
                IOobject
                (
                    "volField",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                 ),
                 mesh,
                 dimensioned<scalar>
                 (
                     "zero",
                     dimVolume,
                     scalar(0.0)
                 ),
                 zeroGradientFvPatchScalarField::typeName
            )
        );
        const scalarField& cvol = mesh.V();

        volScalarField& vols = tVols.ref();
        vols.primitiveFieldRef() = cvol;
        vols.correctBoundaryConditions();

        scalarField neiVols(mesh.nCells(), 0.0);

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();
        scalarField nNei(cvol.size(), 0.0);

        forAll(owner, facei)
        {
            const label own = owner[facei];
            const label nei = neighbour[facei];
            neiVols[own] += cvol[nei];
            neiVols[nei] += cvol[own];
            nNei[own] += 1;
            nNei[nei] += 1;
        }
        forAll(mesh.boundary(), pI)
        {
            if (vols.boundaryField()[pI].coupled())
            {
                tmp<scalarField> tNeiBVols
                (
                    vols.boundaryField()[pI].patchNeighbourField()
                );
                const labelList& pfc = mesh.boundary()[pI].faceCells();
                forAll(pfc, cI)
                {
                    neiVols[pfc[cI]] += tNeiBVols()[cI];
                    nNei[pfc[cI]] += 1;
                }
            }
        }

        forAll(neiVols, cI)
        {
            if (nNei[cI]>0) neiVols[cI] /= nNei[cI];
        }

        scalarField& vol = const_cast<scalarField&>(cvol);

        label nChanged(0);

        forAll(vol, cI)
        {
            scalar& volI = vol[cI];
            const scalar& neivolI = neiVols[cI];
            if (neivolI > 0)
            {
                scalar volFrac = volI/neivolI;
                if (volFrac<threshold)
                {
                    vol[cI] = neivolI*threshold;
                    nChanged++;
                }
            }
            else
            {
                //Special orphan cells: no neis or all zero cells
                if (volI < VSMALL)
                {
                    //- Maybe needs improvement
                    volI = VSMALL;
                }
            }
        }

        if (printInformation)
        {
            reduce(nChanged, sumOp<label>());
            Info<< "Limited cells volume: " << nChanged << endl;
        }
    }
}

// ************************************************************************* //
