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

#include "solidRotationSource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcReconstruct.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "cfdTools/general/findRefCell/findRefCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidRotationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        solidRotationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::solidRotationSource::addLocalFvSchemes
(
    const wordList& fieldNames
)
{
    auto& refMesh = const_cast<fvMesh&>(mesh_);
    dictionary fvSchemesDict = mesh_.schemes().localSchemeDict();
    fvSchemesDict.add("divSchemes", dictionary(), false);
    fvSchemesDict.add("gradSchemes", dictionary(), false);
    fvSchemesDict.add("laplacianSchemes", dictionary(), false);
    for (const word& name : fieldNames)
    {
        fvSchemesDict.subDict("divSchemes").add
        (
            string("div(phiOmega,"+name+")"),
            string("Gauss linearUpwind limitedGrad("+name+")").c_str()
        );
        fvSchemesDict.subDict("gradSchemes").add
        (
            string("limitedGrad("+name+")"),
            "cellLimited Gauss linear 1"
        );
    }
    fvSchemesDict.subDict("laplacianSchemes").add
    (
        string("laplacian(rotatingFluxPotential)"),
        "Gauss linear limited 0.333"
    );
    fvSchemesDict.subDict("gradSchemes").add
    (
        string("grad(rotatingFluxPotential)"),
        "cellLimited Gauss linear 1"
    );
    fvSchemesDict.add("fluxRequired", dictionary(), false);
    fvSchemesDict.subDict("fluxRequired").add("rotatingFluxPotential", word::null);
    refMesh.schemes().setLocalSchemeDict(fvSchemesDict);

    dictionary fvSolutionDict = mesh_.solution().localSolutionDict();
    dictionary solDict
    (
        IStringStream
        (
            "solver GAMG;"
            "smoother symGaussSeidel;"
            "agglomerator faceAreaPair;"
            "mergeLevels 1;"
            "cacheAgglomeration true;"
            "nCellsInCoarsestLevel 20;"
            "scaleCorrection true;"
        )()
    );
    solDict.add("relTol", 1e-3);
    solDict.add("tolerance", 1e-6);
    fvSolutionDict.add("solvers", dictionary(), false);
    fvSolutionDict.subDict("solvers").add
    (
        "rotatingFluxPotential", solDict, false
    );
    solDict.add("relTol", 0, true);
    solDict.add("tolerance", 1e-5, true);
    fvSolutionDict.subDict("solvers").add
    (
        "rotatingFluxPotentialFinal", solDict, false
    );

    refMesh.solution().setLocalSolutionDict(fvSolutionDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidRotationSource::solidRotationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    GRFSource(name, modelType, dict, obr),
    p_
    (
        "rotatingFluxPotential",
        obr,
        mesh(),
        dimVelocity*dimLength,
        0,
        fixedGradientFvPatchScalarField::typeName
    ),
    // Store in object registry as it needs to respond to topo change
    phi0_
    (
        "uncorrectedRotatingFlux",
        obr,
        mesh(),
        dimArea*dimVelocity,
        0
    ),
    nNonOrthCorr_(2),
    writeRotationalVelocity_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::solidRotationSource::read(const dictionary& dict)
{
    if (GRFSource::read(dict))
    {
        coeffs_.readIfPresent("nNonOrthogonalCorrectors", nNonOrthCorr_);
        coeffs_.readIfPresent
        (
            "writeRotationalVelocity", writeRotationalVelocity_
        );
        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::solidRotationSource::phiOmega()
const
{
    tmp<surfaceScalarField> phi0(-GRFSource::phiOmega());
    // First calculation, or the underlying uncorrected flux has changed?
    if (!phi_.valid() || sumMag(phi0()-phi0_).value() > SMALL)
    {
        phi0_ = phi0;
        calculatePhi();
    }

    return tmp<surfaceScalarField>(phi_);
}


void Foam::fv::solidRotationSource::calculatePhi() const
{
    // Correct the rotational flux by disallowing flow through boundaries and
    // imposing a divergence-free condition to interpolate between them

    volScalarField::Boundary& pb = p_.boundaryFieldRef();
    forAll(pb, patchi)
    {
        if (isA<fixedGradientFvPatchScalarField>(pb[patchi]))
        {
            fixedGradientFvPatchScalarField& fvpf =
                refCast<fixedGradientFvPatchScalarField>(pb[patchi]);
            fvpf.gradient() = -phi0_.boundaryField()[patchi];
            if (mesh().moving())
            {
                fvpf.gradient() += mesh().phi().boundaryField()[patchi];
            }
            fvpf.gradient() *=
                this->fvmMask_.boundaryField()[patchi]
               /mesh().magSf().boundaryField()[patchi];
        }
    }

    tmp<fvScalarMatrix> pEqn;
    for (label i = 0; i <= nNonOrthCorr_; i++)
    {
        pEqn = -fvm::laplacian(p_) - this->fvmMask_*fvc::div(phi0_);
        pEqn->setReference(getFirstCell(mesh_).first(), 0);
        pEqn->solve(mesh().solution().solver(p_.select(i == nNonOrthCorr_)));
    }

    phi_.reset((phi0_ - pEqn->flux()).ptr());

    if (writeRotationalVelocity_)
    {
        // It's not necessary to store this, but this is done in order for it to
        // be written together with all the other fields
        rotationalVelocity_.reset
        (
            new volVectorField
            (
                "rotationalVelocity",
                fvc::reconstruct(phi_())
            )
        );
        rotationalVelocity_().writeOpt() = IOobject::AUTO_WRITE;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
