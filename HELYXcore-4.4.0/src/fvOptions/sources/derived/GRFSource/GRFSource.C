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
    (c) 2018-2024 Engys Ltd.
    (c) 2012-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "GRFSource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(GRFSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        GRFSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::GRFSource::addLocalFvSchemes(const wordList& fieldNames)
{
    auto& refMesh = const_cast<fvMesh&>(mesh_);
    dictionary fvSchemesDict = mesh_.schemes().localSchemeDict();
    char scheme[] = ("bounded Gauss limitedLinear 1");
    fvSchemesDict.add("divSchemes", dictionary(), false);
    for (const word& name : fieldNames)
    {
        fvSchemesDict.subDict("divSchemes").add
        (
            string("div(phiOmega,"+name+")"),scheme
        );
        if (name == heTName_)
        {
            fvSchemesDict.subDict("divSchemes").add
            (
                string("div(phiOmega,K)"),scheme
            );
        }
    }

    refMesh.schemes().setLocalSchemeDict(fvSchemesDict);
}


void Foam::fv::GRFSource::setFvmMask()
{
    //set mask field
    const labelListList cellCellLabels = mesh().cellCells();
    forAll(GRFPtrList_, frI)
    {
        forAll(GRFPtrList_[frI]->cells(), cI)
        {
            const label& cL = GRFPtrList_[frI]->cells()[cI];
            fvmMask_[cL] = 1;
            forAll(cellCellLabels[cL], cII)
            {
                const label& cLL = cellCellLabels[cL][cII];
                fvmMask_[cLL] = 1;
            }

        }
    }
    fvmMask_.correctBoundaryConditions();
}


void Foam::fv::GRFSource::setFields() const
{
    phiOmega_ *= 0;
    forAll(GRFPtrList_, frI)
    {
        GRFPtrList_[frI]->applyRotationalFlux(phiOmega_);
    }
}


void Foam::fv::GRFSource::relaxGRFBoundary
(
    fvBlockMatrix<vector>& eqn
) const
{
    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->relaxBoundary(eqn, relax_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::GRFSource::GRFSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(name, modelType, dict, obr),
    GRFPtrList_(0),
    rhoName_(coeffs_.lookupOrDefault<word>("rhoName", "rho")),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    heTName_(coeffs_.lookupOrDefault<word>("heTName", word::null)),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    thermoPtr_(nullptr),
    fvmMask_
    (
        IOobject
        (
            "mrfMask",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    ),
    phiOmega_
    (
        IOobject
        (
            "phiOmega",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol/dimTime, 0)
    ),
    relativePhi_
    (
        IOobject
        (
            "relativePhi",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol/dimTime, 0)
    ),
    consistentDiscretization_
    (
        dict_.lookupOrDefault("consistentDiscretization", false)
    ),
    relax_(dict_.lookupOrDefault<scalar>("relax", 0.5)),
    isInitialised_(false)
{
    MRF_ = true;

    dictionary RFList = coeffs_.subDict("referenceFrames");
    GRFPtrList_.setSize(RFList.toc().size());
    forAll(RFList.toc(), i)
    {
        dictionary RFDict = RFList.subDict(RFList.toc()[i]);
        GRFPtrList_[i].reset(new GRF(RFList.toc()[i], modelType, RFDict, mesh_));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::GRFSource::initialise()
{
    if (!isInitialised_)
    {
        if (obr_.foundObject<volVectorField>(UName_))
        {
            volVectorField& U = obr_.lookupObjectRef<volVectorField>(UName_);
            forAll(U.boundaryField(), patchi)
            {
                if (isA<referenceFrameFvPatch<vector>>(U.boundaryField()[patchi]))
                {
                    const referenceFrameFvPatch<vector>& refPatch =
                        dynamic_cast<const referenceFrameFvPatch<vector>&>
                        (
                            U.boundaryField()[patchi]
                        );
                    refPatch.updateCordinateFrameRegistry();
                    U.boundaryFieldRef()[patchi].resetUpdate();
                }
            }
            U.correctBoundaryConditions();
        }

        thermoPtr_ = obr_.lookupObjectPtr<basicThermo>(basicThermo::dictName);
        if (thermoPtr_)
        {
            heTName_ = thermoPtr_->heT().name();
            TName_ = thermoPtr_->T().name();
        }

        setFvmMask();
        setFields();
        isInitialised_ = true;
    }

    return true;
}


void Foam::fv::GRFSource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames = coeffs_.lookup<wordList>("fields");
    // Replace 'T' with the actual field being solved so the user can
    // just specify 'T' for any energy formulation
    label i = fieldNames.find(TName_);
    if (i >= 0)
    {
        fieldNames[i] = heTName_;
    }
    addLocalFvSchemes(fieldNames);
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::GRFSource::phiOmega() const
{
    setFields();
    phiOmega_.setOriented(true);
    return phiOmega_;
}


void Foam::fv::GRFSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addRotation(geometricOneField(), eqn);
}


void Foam::fv::GRFSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    tmp<surfaceScalarField> rhof;
    if (eqn.psi().name() == "p" || eqn.psi().name() == "p_rgh")
    {
        rhof = fvc::interpolate(rho);
        scalarField& pSource = eqn.source();
        scalarField source(fvc::div(rhof()*phiOmega()));
        const scalarField& V = mesh_.V();
        forAll(mesh().cells(), cI)
        {
            pSource[cI] -= source[cI]*V[cI];
        }
    }
    else if (thermoPtr_ && eqn.psi().name() == TName_)
    {
        // Temperature formulation of energy equation. 'rho' parameter in this
        // case is actually rho*Cv
        const volScalarField& Cv = thermoPtr_->Cv();
        rhof = fvc::interpolate(rho/Cv);
        fvScalarMatrix eqnT(eqn.psi(), dimMass*dimTemperature/dimTime);
        addRotation(rhof(), eqnT);
        eqn += Cv*eqnT;
    }
    else
    {
        rhof = fvc::interpolate(rho);
        addRotation(rhof(), eqn);
    }

    // The energy equation may not be solving total energy, so check for
    // existence of K before applying this correction
    if
    (
        eqn.psi().name() == heTName_
     && eqn.psi().db().foundObject<volScalarField>("K")
    )
    {
        if (obr_.foundObject<surfaceScalarField>("phi"))
        {
            const surfaceScalarField& phi =
                obr_.template lookupObject<surfaceScalarField>("phi");
            const volVectorField& U = obr_.lookupObject<volVectorField>("U");

            surfaceScalarField relativePhi(phi - (rhof*phiOmega()));
            const scalarField& V = mesh_.V();
            scalarField& heTsource = eqn.source();

            scalarField Ksource
            (
                fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")
              - fvc::div(relativePhi, 0.5*magSqr(U), "div(phiOmega,K)")
            );
            if (thermoPtr_ && eqn.psi().name() == TName_)
            {
                Ksource *= thermoPtr_->Cv();
            }
            forAll(mesh().cells(), cI)
            {
                heTsource[cI] -= Ksource[cI]*V[cI];
            }
        }
    }
}


void Foam::fv::GRFSource::addSup
(
    fvVectorMatrix& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(rhoName_);

        forAll(GRFPtrList_, framei)
        {
            GRFPtrList_[framei]->addAcceleration(rho, eqn, true);
        }

        surfaceScalarField rhof(fvc::interpolate(rho));
        addRotation(rhof, eqn);
    }
    else
    {
        forAll(GRFPtrList_, framei)
        {
            GRFPtrList_[framei]->addAcceleration(eqn, true);
        }

        addRotation(geometricOneField(), eqn);
    }
}


void Foam::fv::GRFSource::addSup
(
    const volScalarField& rho,
    fvVectorMatrix& eqn,
    const label fieldi
)
{
    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->addAcceleration(rho, eqn, true);
    }

    surfaceScalarField rhof(fvc::interpolate(rho));
    addRotation(rhof, eqn);
}


void Foam::fv::GRFSource::addAcceleration
(
    fvBlockMatrix<vector>& eqn
) const
{
    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->addAcceleration(eqn);
    }
    addRotation(geometricOneField(), eqn);
}


void Foam::fv::GRFSource::addAcceleration
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn
) const
{
    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->addAcceleration(rho, eqn);
    }
    surfaceScalarField rhof(fvc::interpolate(rho));
    addRotation(rhof, eqn);
}


void Foam::fv::GRFSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    surfaceScalarField rhof(fvc::interpolate(rho));
    addRotation(rhof, eqn);
}


void Foam::fv::GRFSource::makeRelative(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->makeRelative(U);
    }
}


void Foam::fv::GRFSource::zero
(
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    forAll(GRFPtrList_, framei)
    {
        GRFPtrList_[framei]->zero(phi);
    }
}


bool Foam::fv::GRFSource::read(const dictionary& dict)
{
    //use this instead of on-construction initialisation
    if (option::read(dict))
    {
        coeffs_.readIfPresent("rhoName", rhoName_);

        forAll(GRFPtrList_, framei)
        {
            GRFPtrList_[framei]->read();
        }
        initialise();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
