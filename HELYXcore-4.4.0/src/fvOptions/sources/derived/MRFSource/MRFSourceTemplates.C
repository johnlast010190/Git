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
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "MRFSource.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrices.H"
#include "eulerianMultiphaseSystem/eulerianMultiphaseSystem.H"
#include "derivedFvPatchFields/pressureVelocity/pressureVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }


    const labelList& internalFaces = frameSourceFaces_.internalFaces();
    tmp<scalarField> rotField = getRotationalFluxIntField(internalFaces);

    forAll(internalFaces, i)
    {
        const label facei = internalFaces[i];
        phi[facei] -= rho[facei]*rotField()[i];
    }
    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phi
) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField::Boundary& UBf =
        phaseSystemPtr_
      ? phaseSystemPtr_->phases().first().U().boundaryField()
      : obr_.lookupObject<volVectorField>(UName_).boundaryField();
    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];

        // TODO: Unification required the pressure velocity boundary isn't
        // derived from referenceFvPatch, hance it has to be handled
        // separatelly here. Better solution would be to add the
        // referenceFvPatch to directionMixedFvPatch which is base class of
        // pressure velocity bc.
        if
        (
            (
                isA<refFvPatch>(UBf[patchi])
             && dynamic_cast<const refFvPatch&>(UBf[patchi]).inletFlux()
            )
         || (
                isA<pressureVelocityFvPatchVectorField>(UBf[patchi])
             && dynamic_cast<const pressureVelocityFvPatchVectorField&>
                (
                    UBf[patchi]
                ).inletFlux()
           )
        )
        {
            tmp<scalarField> rotBFlux
            (
                getRotationalFluxBField(includedFaces, patchi)
            );
            forAll(includedFaces, i)
            {
                const label patchFacei = includedFaces[i];
                phi[patchi][patchFacei] -=
                    rho[patchi][patchFacei]*rotBFlux()[i];
            }
        }
        else
        {
            UIndirectList<scalar>(phi[patchi], includedFaces) = 0.0;
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        tmp<scalarField> rotBFlux
        (
            getRotationalFluxBField(excludedFaces, patchi)
        );
        forAll(excludedFaces, i)
        {
            label patchFacei = excludedFaces[i];
            phi[patchi][patchFacei] -=
                rho[patchi][patchFacei]*rotBFlux()[i];
        }
    }
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    Field<scalar>& phi,
    const label patchi
) const
{
    if (!active_)
    {
        return;
    }

    const labelList& includedFaces = frameSourceFaces_.includedFaces()[patchi];
    UIndirectList<scalar>(phi, includedFaces) = 0.0;

    const labelList& excludedFaces = frameSourceFaces_.excludedFaces()[patchi];
    forAll(excludedFaces, i)
    {
        const label patchFacei = excludedFaces[i];
        phi[patchFacei] -=
            rho[patchFacei]*getBoundaryRotationalFlux(patchi, patchFacei);
    }
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const labelList& internalFaces = frameSourceFaces_.internalFaces();
    tmp<scalarField> rotField = getRotationalFluxIntField(internalFaces);

    forAll(internalFaces, i)
    {
        const label facei = frameSourceFaces_.internalFaces()[i];
        phi[facei] += rho[facei]*rotField()[i];
    }

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& patchFaces =
            frameSourceFaces_.includedFaces()[patchi];
        tmp<scalarField> rotBFlux
        (
            getRotationalFluxBField(patchFaces, patchi)
        );
        forAll(patchFaces, i)
        {
            const label patchFacei = patchFaces[i];
            scalarField& pphi = phi.boundaryFieldRef()[patchi];
            pphi[patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]*rotBFlux()[i];
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        tmp<scalarField> rotBFlux
        (
            getRotationalFluxBField(excludedFaces, patchi)
        );
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            scalarField& pphi = phi.boundaryFieldRef()[patchi];
            pphi[patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]*rotBFlux()[i];
        }
    }
}


// ************************************************************************* //
