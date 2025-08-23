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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pyrolysisModels/reactingOneDim/reactingOneDim.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvm/fvm.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"
#include "finiteVolume/fvc/fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(reactingOneDim, 0);
    addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, mesh);
    addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::regionModels::reactingOneDim::readReactingOneDimControls()
{
    minimumDelta_ = coeffs().lookup<scalar>("minimumDelta");
    qrHSource_ = coeffs().lookup<bool>("qrHSource");
    useChemistrySolvers_ =
        coeffs().lookupOrDefault<bool>("useChemistrySolvers", true);
}


bool Foam::regionModels::reactingOneDim::read()
{
    if (pyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::regionModels::reactingOneDim::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::regionModels::reactingOneDim::updateqr()
{
    // Update local qr from coupled qr field
    qr_.forceAssign(dimensionedScalar(qr_.dimensions(), 0));

    // Retrieve field from coupled region using mapped boundary conditions
    qr_.correctBoundaryConditions();

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        // qr is positive going in the solid
        // If the surface is emitting the radiative flux is set to zero
        qrBf[patchi] = max(qrBf[patchi], scalar(0));
    }

    const vectorField& cellC = regionMesh().cellCentres();

    tmp<volScalarField> kappa = radiation_().absorptionEmission().a();

    // Propagate qr through 1-D regions
    label localPyrolysisFacei = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        const scalarField& qrp = qr_.boundaryField()[patchi];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchi];

        forAll(qrp, facei)
        {
            const scalar qr0 = qrp[facei];
            point Cf0 = Cf[facei];
            const labelList& cells = boundaryFaceCells_[localPyrolysisFacei++];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label celli = cells[k];
                const point& Cf1 = cellC[celli];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa()[celli]*delta;
                qr_[celli] = qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }
}


void Foam::regionModels::reactingOneDim::topoChange(const scalarField& deltaV)
{
    Info<< "Initial/final volumes = " << gSum(deltaV) << endl;

    // Move the mesh
    const labelList moveMap = moveMesh(deltaV, minimumDelta_);

    // Flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 1)
        {
            solidChemistry_().setCellReacting(i, false);
        }
    }
}


void Foam::regionModels::reactingOneDim::solveContinuity()
{
    if (!moveMesh_)
    {
        fvScalarMatrix rhoEqn(fvm::ddt(rho()) == -solidChemistry_().RRg());
        rhoEqn.solve();
    }
    else
    {
        const scalarField deltaV
        (
            -solidChemistry_().RRg()*regionMesh().V()*time_.deltaT()/rho()
        );

        topoChange(deltaV);
    }
}


void Foam::regionModels::reactingOneDim::solveSpeciesMass()
{
    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn(fvm::ddt(rho(), Yi) == solidChemistry_().RRs(i));

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho())*regionMesh().phi()
            );

            YiEqn -= fvc::div(phiYiRhoMesh);
        }
        YiEqn.solve(regionMesh().solution().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }
    Ys_[Ys_.size() - 1] = 1.0 - Yt;
}


void Foam::regionModels::reactingOneDim::solveEnergy()
{
    volScalarField& he = (*solidThermo_).he();
    tmp<volScalarField> talpha((*solidThermo_).kappa()/(*solidThermo_).Cp());
    volScalarField& alpha = talpha.ref();
    alpha.rename(IOobject::groupName("alpha", he.group()));

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho(), he)
      - fvm::laplacian(alpha, he) + fvc::laplacian(alpha, he)
      - fvc::laplacian((*solidThermo_).kappa(), (*solidThermo_).T())
     ==
        chemistryQdot_
      + solidChemistry_().RRsHs()
    );

    if (qrHSource_)
    {
        const surfaceScalarField phiqr(fvc::interpolate(qr_)*nMagSf());
        heEqn += fvc::div(phiqr);
    }

    heEqn.relax();
    heEqn.solve();
}


void Foam::regionModels::reactingOneDim::calculateMassTransfer()
{
    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistryQdot_);
        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_().RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_().RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::reactingOneDim::reactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    Ys_
    (
        dynamic_cast<solidMulticomponentThermo&>
        (
            *solidThermo_
        ).composition().Y()
    ),
    minimumDelta_(1e-4),
    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime, 0)
    ),
    phiHsGas_
    (
        IOobject("phiHsGas", time().timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime, 0)
    ),
    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    lostSolidMass_(dimensionedScalar(dimMass, 0)),
    addedGasMass_(dimensionedScalar(dimMass, 0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar(dimEnergy/dimTime, 0)),
    qrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active())
    {
        reactingOneDim::read();
    }
}


Foam::regionModels::reactingOneDim::reactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    Ys_
    (
        dynamic_cast<solidMulticomponentThermo&>
        (
            *solidThermo_
        ).composition().Y()
    ),
    minimumDelta_(1e-4),
    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime, 0)
    ),
    phiHsGas_
    (
        IOobject("phiHsGas", time().timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime, 0)
    ),
    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    lostSolidMass_(dimensionedScalar(dimMass, 0)),
    addedGasMass_(dimensionedScalar(dimMass, 0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar(dimEnergy/dimTime, 0)),
    qrHSource_(false),
    useChemistrySolvers_(true)
{
    if (active())
    {
        reactingOneDim::read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::reactingOneDim::~reactingOneDim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::regionModels::reactingOneDim::addMassSources
(
    const label patchi,
    const label facei
)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchi)
        {
            index = i;
            break;
        }
    }

    const label localPatchId = intCoupledPatchIDs_[index];
    const scalar massAdded = phiGas_.boundaryField()[localPatchId][facei];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


Foam::scalar Foam::regionModels::reactingOneDim::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            sqr(regionMesh().surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate((*solidThermo_).kappa())
           /fvc::interpolate((*solidThermo_).Cp()*rho())
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return DiNum;
}


const Foam::surfaceScalarField&
Foam::regionModels::reactingOneDim::phiGas() const
{
    return phiGas_;
}


void Foam::regionModels::reactingOneDim::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    if (useChemistrySolvers_)
    {
        solidChemistry_().solve(time().deltaTValue());
    }
    else
    {
        solidChemistry_().calculate();
    }

    solveContinuity();

    chemistryQdot_ = solidChemistry_().Qdot()();

    if (qrHSource_)
    {
        updateqr();
    }

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=pimple_.correctNonOrthogonal(); nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    (*solidThermo_).correct();

    Info<< "pyrolysis min/max(T) = "
        << min((*solidThermo_).T().primitiveField()) << ", "
        << max((*solidThermo_).T().primitiveField())
        << endl;
}


void Foam::regionModels::reactingOneDim::info()
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << nl
        << indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// ************************************************************************* //
