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

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchFields.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
void Foam::LimitedScheme<Type, Limiter, LimitFunc>::calcLimiter
(
    const VolField<Type>& phi,
    surfaceScalarField& limiterField
) const
{
    typedef VolField<typename Limiter::phiType> VolFieldType;

    typedef VolField<typename Limiter::gradPhiType> GradVolFieldType;

    const fvMesh& mesh = this->mesh();

    tmp<VolFieldType> tlPhi = LimitFunc<Type>()(phi);
    const VolFieldType& lPhi = tlPhi();

    tmp<GradVolFieldType> tgradc(fvc::grad(lPhi));
    const GradVolFieldType& gradc = tgradc();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = limiterField.primitiveFieldRef();

    forAll(pLim, face)
    {
        label own = owner[face];
        label nei = neighbour[face];

        pLim[face] = Limiter::limiter
        (
            CDweights[face],
            this->faceFlux_[face],
            lPhi[own],
            lPhi[nei],
            gradc[own],
            gradc[nei],
            C[nei] - C[own]
        );
    }

    surfaceScalarField::Boundary& bLim = limiterField.boundaryFieldRef();

    forAll(bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            const Field<typename Limiter::phiType> plPhiP
            (
                lPhi.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::phiType> plPhiN
            (
                lPhi.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<typename Limiter::gradPhiType> pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::gradPhiType> pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(CDweights.boundaryField()[patchi].patch().delta());

            forAll(pLim, face)
            {
                pLim[face] = Limiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    plPhiP[face],
                    plPhiN[face],
                    pGradcP[face],
                    pGradcN[face],
                    pd[face]
                );
            }
        }
        else
        {
            pLim = 1.0;
        }
    }

    limiterField.setOriented();
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
Foam::tmp<Foam::surfaceScalarField>
Foam::LimitedScheme<Type, Limiter, LimitFunc>::limiter
(
    const VolField<Type>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    const word limiterFieldName(type() + "Limiter(" + phi.name() + ')');

    if (this->mesh().solution().cache("limiter"))
    {
        if (!mesh.foundObject<surfaceScalarField>(limiterFieldName))
        {
            surfaceScalarField* limiterField
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        limiterFieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimless,
                    fvsPatchField<Type>::calculatedType(),
                    phi.boundaryField().patchTypes()
                )
            );

            mesh.objectRegistry::store(limiterField);
        }

        surfaceScalarField& limiterField =
            const_cast<surfaceScalarField&>
            (
                mesh.lookupObject<surfaceScalarField>(limiterFieldName)
            );

        calcLimiter(phi, limiterField);

        return limiterField;
    }
    else
    {
        tmp<surfaceScalarField> tlimiterField
        (
            surfaceScalarField::New
            (
                limiterFieldName,
                mesh,
                dimless,
                fvsPatchField<Type>::calculatedType(),
                phi.boundaryField().patchTypes()
            )
        );

        calcLimiter(phi, tlimiterField.ref());

        return tlimiterField;
    }
}


// ************************************************************************* //
