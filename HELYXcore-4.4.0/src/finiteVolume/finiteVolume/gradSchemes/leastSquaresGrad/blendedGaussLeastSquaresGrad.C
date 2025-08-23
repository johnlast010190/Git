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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/leastSquaresGrad/blendedGaussLeastSquaresGrad.H"
#include "finiteVolume/gradSchemes/leastSquaresGrad/leastSquaresVectors.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::blendedGaussLeastSquaresGrad<Type>::calcGrad
(
    const VolField<Type>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<VolField<GradType>> tlsGrad
    (
        VolField<GradType>::New
        (
            name,
            vsf.db(),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    VolField<GradType>& lsGrad = tlsGrad.ref();

    // Get reference to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();
    const volScalarField layerCells = lsv.wallLayerCells(nWallLaynerCells_);

    const surfaceVectorField& Sf = mesh.Sf();
    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
    const surfaceVectorField& mSf = tmSf();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    tmp<SurfaceField<Type>> tssf
    (
        linearInterpolate(vsf)
    );
    const SurfaceField<Type>& issf = tssf();

    forAll(own, facei)
    {
        label ownFacei = own[facei];
        label neiFacei = nei[facei];

        Type deltaVsf = vsf[neiFacei] - vsf[ownFacei];
        GradType Sfssf = mSf[facei]*issf[facei];
        if (layerCells[ownFacei] == 1)
        {
            lsGrad[ownFacei] += Sfssf;
        }
        else
        {
            lsGrad[ownFacei] += ownLs[facei]*deltaVsf;
        }
        if (layerCells[neiFacei] == 1)
        {
            lsGrad[neiFacei] -= Sfssf;
        }
        else
        {
            lsGrad[neiFacei] -= neiLs[facei]*deltaVsf;
        }
    }

    // Boundary faces

    // treatment for extrapolated boundaries
    autoPtr<tensorField> gradScalePtr(nullptr);
    const scalarField& Vol(mesh.V());

    forAll(vsf.boundaryField(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);

        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const labelUList& faceCells =
            vsf.boundaryField()[patchi].patch().faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        if (vsf.boundaryField()[patchi].extrapolated())
        {
            if (!gradScalePtr.valid())
            {
                gradScalePtr.reset(new tensorField(vsf.size(), I));
            }

            tmp<vectorField> delta(patch.Cf() - patch.Cn());

            forAll(patch, patchFaceI)
            {
                const label celli = faceCells[patchFaceI];
                if (layerCells[celli] == 1)
                {
                    //scaling contribution
                    gradScalePtr->operator[](celli)
                        -= (pmSf[patchFaceI]*delta->operator[](patchFaceI))/Vol[celli];

                    //zero gradient contributions
                    lsGrad[celli] += (pmSf[patchFaceI]*vsf[celli]);
                }
                else
                {
                    gradScalePtr->operator[](celli) -=
                        (
                            patchOwnLs[patchFaceI]
                           *delta->operator[](patchFaceI)
                        );
                }
            }
        }
        else if (vsf.boundaryField()[patchi].coupled())
        {
            const tmp<Field<Type>> neiVsf
            (
                vsf.boundaryField()[patchi].patchNeighbourField()
            );
            const fvsPatchField<Type>& pssf = tssf->boundaryField()[patchi];

            forAll(neiVsf(), patchFacei)
            {
                const label celli = faceCells[patchFacei];
                if (layerCells[celli] == 1)
                {
                    lsGrad[celli] += pmSf[patchFacei]*pssf[patchFacei];
                }
                else
                {
                    lsGrad[celli] +=
                        patchOwnLs[patchFacei]
                       *(neiVsf()[patchFacei] - vsf[celli]);
                }
            }
        }
        else
        {
            tmp<Field<Type>> tpatchVsf =
                vsf.boundaryField()[patchi].gradientBoundaryValue();
            const Field<Type>& patchVsf = tpatchVsf();

            forAll(patchVsf, patchFacei)
            {
                const label celli = faceCells[patchFacei];
                if (layerCells[celli] == 1)
                {
                    lsGrad[celli]
                        += pmSf[patchFacei]*patchVsf.operator[](patchFacei);
                }
                else
                {
                    lsGrad[celli] +=
                         patchOwnLs[patchFacei]
                        *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
                }
            }
        }
    }

    //since multiple extrapolation boundaries can contribute to the same cell
    //the multiplication has top be done after the summation is complete
    if (gradScalePtr.valid())
    {
        forAll(gradScalePtr(), celli)
        {
            if (gradScalePtr->operator[](celli) != tensor::I)
            {
                lsGrad[celli]
                    = (inv(gradScalePtr->operator[](celli)) & lsGrad[celli]);
            }
        }
    }
    forAll(layerCells, celli)
    {
        if (layerCells[celli] == 1)
        {
            lsGrad[celli] /= Vol[celli];
        }
    }


    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


template<class Type>
Foam::tmp
<
    Foam::BlockLduSystem
    <
        Foam::vector, typename Foam::outerProduct<Foam::vector, Type>::type
    >
> Foam::fv::blendedGaussLeastSquaresGrad<Type>::fvmGrad
(
    const VolField<Type>& vf
) const
{
    FatalErrorInFunction
        << "Implicit gradient operator defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType>> tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// ************************************************************************* //
