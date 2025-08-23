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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "MRFSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "eulerianMultiphaseSystem/eulerianMultiphaseSystem.H"
#include "fvMesh/fvPatches/constraint/nonConformal/nonConformalFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(MRFSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        MRFSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalar Foam::fv::MRFSource::getRotationalFlux
(
    const label facei
) const
{
    if (facei >= mesh_.nInternalFaces())
    {
        FatalErrorInFunction
            << "Boundary face treated as internal face."
            << exit(FatalError);
    }

    const face& f = mesh_.faces()[facei];
    if (f.size() > 3 && conservative())
    {
        return getTriangularFlux(f);
    }
    else
    {
        return
            coorFrame().frameVelocity(mesh_.Cf()[facei], false)
          & mesh_.Sf()[facei];
    }
}


Foam::tmp<Foam::scalarField> Foam::fv::MRFSource::getRotationalFluxIntField
(
    const labelList& internalFaces
) const
{
    tmp<scalarField> tRotFlux(new scalarField(internalFaces.size(), Zero));

    if (!conservative())
    {
        vectorField cfInt(mesh_.Cf()(), internalFaces);
        vectorField sfInt(mesh_.Sf()(), internalFaces);
        tRotFlux.ref() = coorFrame().frameVelocity(cfInt, false) & sfInt;
    }
    else
    {
        DynamicList<label>  dmap(4*internalFaces.size());
        DynamicList<vector> dcf(4*internalFaces.size());
        DynamicList<vector> dsf(4*internalFaces.size());
        const pointField& p = mesh_.points();
        forAll(internalFaces, fI)
        {
            const label facei = internalFaces[fI];
            const face& f = mesh_.faces()[facei];
            if (f.size() > 3)
            {
                label nPoints = f.size();
                point fCentre = p[f[0]];
                for (label pi = 1; pi < nPoints; pi++)
                {
                    fCentre += p[f[pi]];
                }
                for (label pi = 0; pi < nPoints; pi++)
                {
                    const point& nextPoint = p[f[(pi + 1) % nPoints]];
                    vector c = p[f[pi]] + nextPoint + fCentre;
                    c /= 3.;
                    vector Sf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                    dmap.append(facei);
                    dcf.append(c);
                    dsf.append(Sf);
                }
            }
            else
            {
                dmap.append(facei);
                dcf.append(mesh_.Cf()[facei]);
                dsf.append(mesh_.Sf()[facei]);
            }
        }
        dmap.shrink();
        dcf.shrink();
        dsf.shrink();
        vectorField cf(dcf);
        vectorField sf(dsf);
        tmp<scalarField> trotFluxCons = coorFrame().frameVelocity(cf, false)&sf;

        forAll(dmap, fI)
        {
            const label fII = dmap[fI];
            tRotFlux.ref()[fII] += trotFluxCons()[fI];
        }
    }

    return tRotFlux;
}


Foam::tmp<Foam::scalarField> Foam::fv::MRFSource::getRotationalFluxBField
(
    const labelList& pFaces,
    const label patchI
) const
{
    tmp<scalarField> tRotFlux(new scalarField(pFaces.size(), Zero));

    bool isConformal = isA<nonConformalFvPatch>(mesh_.boundary()[patchI]);

    if (!conservative() || isConformal)
    {
        vectorField cfb(mesh_.Cf().boundaryField()[patchI], pFaces);
        vectorField sfb(mesh_.Sf().boundaryField()[patchI], pFaces);
        tRotFlux.ref() = coorFrame().frameVelocity(cfb, false) & sfb;
    }
    else
    {
        DynamicList<label>  dmap(4*pFaces.size());
        DynamicList<vector> dcf(4*pFaces.size());
        DynamicList<vector> dsf(4*pFaces.size());
        const pointField& p = mesh_.points();
        label startFace = mesh_.boundaryMesh()[patchI].start();
        forAll(pFaces, fI)
        {
            const label facei = pFaces[fI];
            const label gfacei = facei + startFace;
            const face& f = mesh_.faces()[gfacei];
            if (f.size() > 3)
            {
                label nPoints = f.size();
                point fCentre = p[f[0]];
                for (label pi = 1; pi < nPoints; pi++)
                {
                    fCentre += p[f[pi]];
                }
                for (label pi = 0; pi < nPoints; pi++)
                {
                    const point& nextPoint = p[f[(pi + 1) % nPoints]];
                    vector c = p[f[pi]] + nextPoint + fCentre;
                    c /= 3.;
                    vector Sf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                    dmap.append(facei);
                    dcf.append(c);
                    dsf.append(Sf);
                }
            }
            else
            {
                dmap.append(facei);
                dcf.append(mesh_.Cf().boundaryField()[patchI][facei]);
                dsf.append(mesh_.Sf().boundaryField()[patchI][facei]);
            }
        }
        dmap.shrink();
        dcf.shrink();
        dsf.shrink();
        vectorField cf(dcf);
        vectorField sf(dsf);
        tmp<scalarField> trotFluxCons = coorFrame().frameVelocity(cf, false)&sf;

        forAll(dmap, fI)
        {
            const label fII = dmap[fI];
            tRotFlux.ref()[fII] += trotFluxCons()[fI];
        }
    }

    return tRotFlux;
}


Foam::scalar Foam::fv::MRFSource::getBoundaryRotationalFlux
(
    const label patchi,
    const label patchFacei
) const
{
    // Only default handling for non-conformal patch faces
    if (isA<nonConformalFvPatch>(mesh_.boundary()[patchi]))
    {
        const vector& Cf = mesh_.Cf().boundaryField()[patchi][patchFacei];
        const vector& Sf = mesh_.Sf().boundaryField()[patchi][patchFacei];

        return (coorFrame().frameVelocity(Cf, false)) & Sf;
    }

    const label facei = patchFacei + mesh_.boundaryMesh()[patchi].start();
    const face& f = mesh_.faces()[facei];
    if (f.size() > 3 && conservative())
    {
        return getTriangularFlux(f);
    }
    else
    {
        const vector& Cf = mesh_.Cf().boundaryField()[patchi][patchFacei];
        const vector& Sf = mesh_.Sf().boundaryField()[patchi][patchFacei];

        return (coorFrame().frameVelocity(Cf, false)) & Sf;
    }
}


Foam::scalar Foam::fv::MRFSource::getTriangularFlux
(
    const face& f
) const
{
    const pointField& p = mesh_.points();
    scalar totalRotFlux = 0;
    point fCentre = p[f[0]];
    label nPoints = f.size();
    for (label pi = 1; pi < nPoints; pi++)
    {
        fCentre += p[f[pi]];
    }
    for (label pi = 0; pi < nPoints; pi++)
    {
        const point& nextPoint = p[f[(pi + 1) % nPoints]];
        vector c = p[f[pi]] + nextPoint + fCentre;
        c /= 3.;
        vector Sf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
        totalRotFlux += (coorFrame().frameVelocity(c, false) & Sf);
    }
    return totalRotFlux;
}


void Foam::fv::MRFSource::updateBoundaryVelocity(volVectorField& U)
{
    forAll(U.boundaryField(), patchi)
    {
        if (isA<refFvPatch>(U.boundaryField()[patchi]))
        {
            const refFvPatch& frameFvPatch =
                dynamic_cast<const refFvPatch&>(U.boundaryField()[patchi]);
            frameFvPatch.updateCordinateFrameRegistry();
            U.boundaryFieldRef()[patchi].resetUpdate();
        }
    }
    U.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MRFSource::MRFSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    frameSources(name, modelType, dict, obr),
    rhoName_(coeffs_.lookupOrDefault<word>("rhoName", "rho")),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    isInitialised_(false),
    phaseSystemPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::MRFSource::initialise()
{
    if (!isInitialised_)
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, coeffs_);
        coorFrame().registry().attachObject(this->name());

        phaseSystemPtr_ =
            obr_.lookupObjectRefPtr<eulerianMultiphaseSystem>
            (
                eulerianPhaseSystem::typeName
            );

        if (phaseSystemPtr_)
        {
            // Do all phasic velocities
            forAllIters(phaseSystemPtr_->phases(), iter)
            {
                updateBoundaryVelocity(iter().U());
            }
        }
        else
        {
            volVectorField& U = obr_.lookupObjectRef<volVectorField>(UName_);
            updateBoundaryVelocity(U);
        }

        MRF_ = true;
        isInitialised_ = true;
    }

    return true;
}


void Foam::fv::MRFSource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames = wordList();
}


void Foam::fv::MRFSource::addSup
(
    fvVectorMatrix& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        addAcceleration
        (
            obr_.lookupObject<volScalarField>(rhoName_),
            eqn,
            false
        );
    }
    else
    {
        addAcceleration(eqn, false);
    }
}


void Foam::fv::MRFSource::addSup
(
    const volScalarField& rho,
    fvVectorMatrix& eqn,
    const label fieldI
)
{
    addAcceleration(rho, eqn, false);
}


void Foam::fv::MRFSource::makeRelative(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        U[celli] -= coorFrame().frameVelocity(C[celli], false);
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();
    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];
        forAll(includedFaces, i)
        {
            const label patchFacei = includedFaces[i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }
}


void Foam::fv::MRFSource::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::fv::MRFSource::makeRelative
(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::fv::MRFSource::makeRelative
(
    Field<scalar>& phi,
    const label patchi
) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::fv::MRFSource::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::fv::MRFSource::makeAbsolute(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        U[celli] += coorFrame().frameVelocity(C[celli], false);
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];
        forAll(includedFaces, i)
        {
            const label patchFacei = includedFaces[i];
            Ubf[patchi][patchFacei] +=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            Ubf[patchi][patchFacei] +=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }
}


void Foam::fv::MRFSource::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::fv::MRFSource::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::fv::MRFSource::zero
(
    surfaceScalarField& phi
) const
{
    if (!active_ || ddtPhiCorr())
    {
        return;
    }

    UIndirectList<scalar>(phi, frameSourceFaces_.internalFaces()) = 0.0;

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& patchFaces =
            frameSourceFaces_.includedFaces()[patchi];
        const scalarField& pphi = phi.boundaryFieldRef()[patchi];
        UIndirectList<scalar>(pphi, patchFaces) = 0.0;
    }
    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        const scalarField& pphi = phi.boundaryFieldRef()[patchi];
        UIndirectList<scalar>(pphi, excludedFaces) = 0.0;
    }
}


void Foam::fv::MRFSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::MRFSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("rhoName", rhoName_);

        initialise();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
