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
    (c) 2011 Chenxiaoxiao
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "driftSedimentationSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "rhoThermo/rhoThermo.H"
#include "fields/fvPatchFields/derived/sedimentation/sedimentationFvPatchField.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "pointPatchFields/derived/sedimentation/sedimentationPointPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(driftSedimentationSolver, 0);
    }
}

makeFvSolverOption(driftSedimentationSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::driftSedimentationSolver::driftSedimentationSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    UResidual_({1, 1}),
    pResidual_({1, 1}),
    fieldResidual_({1, 1}),
    flowSolutionFinished_(false),
    fieldSolved_(false),
    turbulencePtr_(nullptr),
    fieldCorrectorFinished_(false)
{
    const Time& runTime = mesh().time();
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    Info<< "Reading thermophysical properties\n" << endl;
    thermoPtr_ =
        &refCast<rhoThermo>
        (
            multiphaseThermo::lookupOrCreate(obr_, phaseName_)
        );

    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    fieldName_ = dict.lookupOrDefault<word>("fieldName", "sediment");

    Info<< "Reading sediment field\n" << endl;
    field_.set
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,
                runTime.timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Info<< "Creating sediment velocity field\n" << endl;
    wf_.set
    (
        new volVectorField
        (
            IOobject
            (
                "Wf",
                runTime.timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("wf", dimVelocity, dict_)
        )
    );

    // Get ids of the sediment patches
    wordList sedimentPatchNames;
    sedimentPatchIDs_.resize(0);
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<sedimentationFvPatchField>(field_().boundaryField()[patchi]))
        {
            sedimentPatchIDs_.append(patchi);
            sedimentPatchNames.append(mesh_.boundaryMesh()[patchi].name());
        }
    }

    // Patches to be included for U internal contribution
    const wordList UInternalPatches =
        dict_.lookupOrDefault<wordList>("UInternalPatches", wordList());
    UInternalPatches_ = boolList(sedimentPatchIDs_.size(), false);
    forAll(sedimentPatchIDs_, i)
    {
        const word& name = sedimentPatchNames[i];
        if (UInternalPatches.found(name))
        {
            UInternalPatches_[i] = true;
        }
    }

    volScalarField& rho =
        obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));
    phiWf_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiWf",
                runTime.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(rho*wf_()) & mesh_.Sf()
        )
    );
    updatePhiWf(true);

    // Unlike the original the field should use READ_IF_PRESENT
    // to ensure the consitent restart of the simulation.
    Info<< "Creating total mass exchange field\n" << endl;
    M_.set
    (
        new volScalarField
        (
            IOobject
            (
                "M",
                runTime.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/dimTime, 0)
        )
    );

    // Unlike the original the field should use READ_IF_PRESENT
    // to ensure the consitent restart of the simulation.
    Info<< "Creating sediment height field\n" << endl;
    deltaH_.set
    (
        new volScalarField
        (
            IOobject
            (
                "deltaH",
                runTime.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, 0)
        )
    );
    deltaHold_.set(new volScalarField("deltaHold", deltaH_));

    Info<< "Creating sediment shear stress field\n" << endl;
    shearStress_.set
    (
        new volVectorField
        (
            IOobject
            (
                "shearStress",
                runTime.timeName(),
                obr,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "shearStress",
                dimMass/dimLength/sqr(dimTime),
                Zero
            )
        )
    );

    // Filds for post-processing
    Info<< "Creating filed with initial cell/face center locations" << endl;

    normal_ =
        -normalised
        (
            uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    mesh().time().constant(),
                    obr_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ).value()
        );

    startingSurfaceName_ =
        dict.lookupOrDefault<word>("startingSurface", word::null);

    // Read or create face flux
    if
    (
        IOobject
        (
            "surfaceLevel",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            false
        ).typeHeaderOk<volVectorField>()
    )
    {
        surfaceLevel_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "surfaceLevel",
                    runTime.timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
    else
    {
        surfaceLevel_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "surfaceLevel",
                    runTime.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("surfaceLevel", dimLength, Zero),
                fixedValueFvPatchField<scalar>::typeName
            )
        );
        autoPtr<triSurfaceMesh> startingSurface;
        if (startingSurfaceName_ != word::null)
        {
            startingSurface.set
            (
                new triSurfaceMesh
                (
                    IOobject
                    (
                        startingSurfaceName_,
                        runTime.constant(),
                        "triSurface",
                        obr_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
        }
        forAll(sedimentPatchIDs_, i)
        {
            const label patchi = sedimentPatchIDs_[i];
            const vectorField& Cfs = mesh_.Cf().boundaryField()[patchi];
            autoPtr<vectorField> nearestPoints;
            if (startingSurface.valid())
            {
                // Find all intersections with the surface
                nearestPoints.set
                (
                    new vectorField
                    (
                        startingSurface().findRayNearest(Cfs, -normal_)
                    )
                );
            }
            vectorField& values =
                surfaceLevel_().boundaryFieldRef()[patchi];
            forAll(values, facei)
            {
                if
                (
                    nearestPoints.valid()
                 && nearestPoints()[facei] != vector(GREAT, GREAT, GREAT)
                 && (
                        (nearestPoints()[facei] & normal_)
                      < (Cfs[facei] & normal_)
                    )
                )
                {
                    values[facei] = nearestPoints()[facei];
                }
                else
                {
                    values[facei] = Cfs[facei];
                }
            }
        }
    }

    // Distance for the U internal
    UInternalDistance_ =
        dict_.lookupOrDefault<scalar>("UInternalDistance", 0.1);

    // TODO: Selection of the cell labels should be properly synced however that
    // requires quite a bit additional code and will be considered only if
    // needed in future. For now the algorithm will use the nearest cell
    // to the patch when cell isn't found at specified distance from the patch.
    forAll(sedimentPatchIDs_, i)
    {
        if (UInternalPatches_[i])
        {
            const label patchi = sedimentPatchIDs_[i];
            const vectorField& Cfs = mesh_.Cf().boundaryField()[patchi];
            const labelUList& cellsNextToPatch =
                shearStress_().boundaryField()[patchi].patch().faceCells();
            labelList cellIds(Cfs.size(), -1);
            forAll(cellIds, j)
            {
                vector offsetPoint = Cfs[j] + normal_*UInternalDistance_;
                cellIds[j] = mesh_.findNearestCell(offsetPoint);
                // If cell isn't found it algorithm will use cell next to the
                // patch
                cellIds[j] =
                    (cellIds[j] != -1) ? cellIds[j] : cellsNextToPatch[j];
            }
            patchesCells_.set
            (
                patchi,
                autoPtr<labelList>(new labelList(cellIds))
            );
        }
    }

    // Total sediment depth will be calculated on "sediment patches"
    totalSedimentDepth_.set
    (
        new volScalarField
        (
            IOobject
            (
                "totalSedimentDepth",
                runTime.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, 0)
        )
    );
    if
    (
        !IOobject
        (
            "totalSedimentDepth",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            false
        ).typeHeaderOk<volVectorField>()
    )
    {
        Info<< "Initialising total sediment depth\n" << endl;
        forAll(sedimentPatchIDs_, i)
        {
            const label patchi = sedimentPatchIDs_[i];
            totalSedimentDepth_().boundaryFieldRef()[patchi] =
                (
                   normal_
                 & (
                       mesh_.Cf().boundaryField()[patchi]
                     - surfaceLevel_().boundaryField()[patchi]
                   )
                );
        }
    }

    sedimentLoad_.set
    (
        new volScalarField
        (
            IOobject
            (
                "sedimentLoad",
                runTime.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/dimArea, 0)
        )
    );

    // sediment density
    rhoSediment_ = dict_.lookupOrDefault<scalar>("rhoSediment", 150);

    // Controlling constant of erosion equation
    erosionConstant_ = dict_.lookupOrDefault<scalar>("erosionConstant", 7e-4);

    // Threshold shear wind speed
    thresholdShearWindSpeed_ =
        dict_.lookupOrDefault<scalar>("Uthreshold", 0.2);

    const dictionary& resDict = dict_.subDict("residuals");

    // Velocity residual threshold in sub-cycles
    UResidual_.input = resDict.lookupOrDefault<scalar>("U", 5e-4);

    // Pressure residual threshold in sub-cycles
    pResidual_.input = resDict.lookupOrDefault<scalar>("p", 5e-4);

    // Concentration residual threshold in sub-cycles
    fieldResidual_.input =
        resDict.lookupOrDefault<scalar>(field_().name(), 5e-6);

    // Maximum sub-cycles per stage
    maxOuterCorrectors_ = dict_.lookupOrDefault<label>("nMaxSubCycles", 1000);

    // Maximum sub-cycles for sediment.
    maxSedimentCorrectors_ =
        dict_.lookupOrDefault<label>("nMaxSedimentSubCycles", 1000);

    alternativeErosion_ =
        dict_.lookupOrDefault<Switch>("alternativeErosion", false);

    // Using tominaga Deposition Rate Phip*Wf*A
    TominagaDepositionRate_ =
        dict_.lookupOrDefault<Switch>("TominagaDepositionRate", false);

    // Using tominaga Deposition Rate Phip*Wf*A
    depositionWithErosion_ =
        dict_.lookupOrDefault<Switch>("depostionWithErosion", false);

    // Report information about loaded variables
    Info<< nl << "========= Setup for the "
        << driftSedimentationSolver::typeName << " =========" << nl
        << "    Treshold velocity for errosion: "
        << thresholdShearWindSpeed_ << nl
        << "    Erosion constant: " << erosionConstant_ << nl
        << "    Sediment density: " << rhoSediment_ << nl
        << "    Sediment patches: " << sedimentPatchNames << nl << nl
        << "Solution controls:" << nl
        << "    nMaxSubCycles: " << maxOuterCorrectors_ << nl
        << "    nMaxSedimentSubCycles: " << maxSedimentCorrectors_ << nl
        << "    residuals" << nl
        << "    {" << nl
        << "        U " << UResidual_.input << nl
        << "        p " << pResidual_.input << nl
        << "        "  << field_().name() << " "
        << fieldResidual_.input << nl
        << "    }"
        << nl << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::driftSedimentationSolver::updatePhiWf(bool updateOnlyBoundary)
{
    volScalarField& rho =
        obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));

    if (!updateOnlyBoundary)
    {
        phiWf_() = linearInterpolate(rho*wf_()) & mesh_.Sf();
    }

    forAll(sedimentPatchIDs_, i)
    {
        const label patchi = sedimentPatchIDs_[i];

        // The falling velocity at the sediment surface is 0.
        phiWf_().boundaryFieldRef()[patchi] = Zero;
    }
}


bool Foam::fv::driftSedimentationSolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    // Turbulence model is compulsory
    turbulencePtr_ =
        &obr_.lookupObject<compressible::turbulenceModel>
        (
            compressible::turbulenceModel::propertiesName
        );

    return true;
}


void Foam::fv::driftSedimentationSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solve always after turbulence
    solveNames = {field_->name()};
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
    optionalDependencies.insert(field_->name(), {"turbulence"});

    // Field correctors
    correctorMembers.insert("corrector:" + field_->name(), {field_->name()});
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::fv::driftSedimentationSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    if
    (
        fieldName == field_->name()
     && flowSolutionFinished_
     && !fieldCorrectorFinished_
    )
    {
        volScalarField& rho =
            obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));

        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>("phi");

        updatePhiWf();

        // Diffusion coefficient
        volScalarField D("D" + field_().name(), turbulencePtr_->Dmt());

        tmp<fvScalarMatrix> tEqn
        (
            new fvScalarMatrix
            (
                fvm::div(phi, field_())
              + fvm::div(phiWf_(), field_())
              - fvm::laplacian(D, field_())
            ==
                fvOptions()(rho, field_())
            )
        );

        tEqn.ref().relax();
        fvOptions().constrain(tEqn.ref());
        return tEqn;
    }

    return tmp<fvScalarMatrix>();
}


void Foam::fv::driftSedimentationSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    if
    (
        solveName == field_->name()
     && flowSolutionFinished_
     && !fieldCorrectorFinished_
    )
    {
        fieldSolved_ = true;
        this->fvOptions().correct(field_());

        // Shear stress/rho in [kg/(m.s^2)]
        tmp<volSymmTensorField> tReff(turbulencePtr_->devRhoReff());
        volSymmTensorField& Reff = tReff.ref();
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(addPhaseName("rho"));

        // Correct the mass exchange rate M.
        forAll(sedimentPatchIDs_, i)
        {
            const label patchi = sedimentPatchIDs_[i];

            // Area vector
            const vectorField& Sfp =
                this->mesh().Sf().boundaryField()[patchi];

            // Magnitude of the area vector
            const scalarField& magSfp =
                this->mesh().magSf().boundaryField()[patchi];

            const scalarField& rhop = rho.boundaryField()[patchi];

            // Concentration of sediment drift at the boundary
            const scalarField& sedimentp =
                field_().boundaryField()[patchi];
            const symmTensorField& Reffp = Reff.boundaryField()[patchi];

            // Shear stress
            vectorField& shearStress =
                shearStress_().boundaryFieldRef()[patchi];
            shearStress = (-Sfp/magSfp) & Reffp;

            // Magnitude of shear velocity [m/s]
            scalarField UShear(sqrt(mag(shearStress/rhop)));
            scalarField& Mp = M_().boundaryFieldRef()[patchi];
            scalarField& deltaHp = deltaH_().boundaryFieldRef()[patchi];
            vectorField& wfp = wf_().boundaryFieldRef()[patchi];
            const volVectorField& U =
                obr_.lookupObjectRef<volVectorField>("U");

            // Update the mass exchange rate (erosion/deposition)
            // update mass exchange rate on sediment surface
            // (erosion & deposition)
            forAll(Mp, facei)
            {
                // A = Ahol => area projected the sedimentation plane
                const scalar area = mag(-Sfp[facei] & normal_);

                // Deposition contribution from internal field
                scalar Udepos = 0.0;
                if (UInternalPatches_[i])
                {
                    const label celli = patchesCells_[patchi]()[facei];
                    scalar UInternalShear =
                        mag(U[celli] - normal_*(normal_ & U[celli]));
                    Udepos = U[celli] & normal_;
                    Udepos =
                        (
                            Udepos > 0.0
                         || UInternalShear > thresholdShearWindSpeed_
                        ) ? 0.0 : Udepos;
                }

                // Reset the value since MpTotal can be Merosion + Mdeposition
                Mp[facei] = 0.0;

                // UShear[facei] > Uthreshold  => Erosion otherwise Deposition
                if (UShear[facei] > thresholdShearWindSpeed_)
                {
                    if (alternativeErosion_)
                    {
                        Mp[facei] =
                            UShear[facei] < SMALL
                          ? 0.0
                          : -erosionConstant_*rhoSediment_*UShear[facei]
                           *(
                                1.0
                              - sqr(thresholdShearWindSpeed_)
                               /sqr(UShear[facei])
                            )
                           *area;
                    }
                    else
                    {
                        Mp[facei] =
                           -erosionConstant_
                           *(sqr(UShear[facei]) - sqr(thresholdShearWindSpeed_))
                           *area;
                    }
                }

                if
                (
                    depositionWithErosion_
                 || UShear[facei] <= thresholdShearWindSpeed_
                )
                {
                    const scalar MpNoShear =
                       -sedimentp[facei]
                       *((wfp[facei] & normal_) + Udepos)*area;

                    if (TominagaDepositionRate_)
                    {
                        Mp[facei] += MpNoShear;
                    }
                    else
                    {
                        Mp[facei] +=
                            MpNoShear
                           *(
                                1.0
                              - sqr(UShear[facei]/thresholdShearWindSpeed_)
                            );
                    }
                }

                deltaHp[facei] =
                    (Mp[facei]/rhoSediment_/area)*mesh().time().deltaTValue();
            }
        }
    }
}


void Foam::fv::driftSedimentationSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::timeLoopCorrectorName)
    {
        const_cast<Time&>(mesh().time()).endSubCycle();

        // Calculate total sediment depth
        forAll(sedimentPatchIDs_, i)
        {
            const label patchi = sedimentPatchIDs_[i];
            scalarField& ptotSedDepth =
                totalSedimentDepth_().boundaryFieldRef()[patchi];

            const pointVectorField* pointMotionU =
                obr_.lookupObjectPtr<pointVectorField>("pointMotionU");

            if
            (
                !pointMotionU
             || !isA<sedimentationPointPatchVectorField>
                (
                    pointMotionU->boundaryField()[patchi]
                )
            )
            {
                ptotSedDepth += deltaHold_().boundaryField()[patchi];

                // Bound the total sediment depth to not go bellow specified
                // geometry constrain
                forAll(ptotSedDepth, facei)
                {
                    ptotSedDepth[facei] = max(ptotSedDepth[facei], 0.0);
                }
            }
            else
            {
                ptotSedDepth =
                    (
                        normal_
                      & (
                            mesh_.Cf().boundaryField()[patchi]
                          - surfaceLevel_().boundaryField()[patchi]
                        )
                    );
            }
        }
        deltaHold_.reset(new volScalarField("deltaHold", deltaH_));

        // Total sediment mass per patch
        Info<< nl;
        forAll(sedimentPatchIDs_, i)
        {
            const label patchi = sedimentPatchIDs_[i];
            sedimentLoad_().boundaryFieldRef()[patchi] =
                totalSedimentDepth_().boundaryField()[patchi]*rhoSediment_;

            scalarField massOnFaces
            (
                sedimentLoad_().boundaryField()[patchi]
               *mesh_.magSf().boundaryField()[patchi]
            );

            Info<< "Total sediment mass on patch \""
                << mesh_.boundaryMesh()[patchi].name() << "\" is: "
                << gSum(massOnFaces) << " kg. " << nl;
        }
    }
}


bool Foam::fv::driftSedimentationSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    Time& runTime = const_cast<Time&>(mesh().time());
    const word& outerCorrName = solverObject::outerCorrectorName;
    if (corrector == 0 && correctorName == outerCorrName)
    {
        // Subcycles needed to update residuals otherwise initialResidual is
        // kept constant within one timestep.
        const scalar deltaT = runTime.deltaTValue();
        runTime.subCycle(maxOuterCorrectors_ + maxSedimentCorrectors_);
        runTime.setDeltaT(deltaT);
        flowSolutionFinished_ = false;
        fieldSolved_ = false;
    }
    else if
    (
        (correctorName == outerCorrName || flowSolutionFinished_)
     && correctorName != solverObject::timeLoopCorrectorName
    )
    {
        runTime++;
    }

    const Foam::dictionary& solverDict = mesh().solverPerformanceDict();

    UResidual_.solution = 1;
    pResidual_.solution = 1;
    fieldResidual_.solution = 1;
    if (solverDict.found("U"))
    {
        UResidual_.solution =
            mag
            (
                List<SolverPerformance<vector>>
                (
                    solverDict.lookup("U")
                ).first().initialResidual()
            );
        pResidual_.solution =
            List<SolverPerformance<scalar>>
            (
                solverDict.lookup("p")
            ).first().initialResidual();

        if (solverDict.found(field_->name()) && fieldSolved_)
        {
            fieldResidual_.solution =
                List<SolverPerformance<scalar>>
                (
                    solverDict.lookup(field_->name())
                ).first().initialResidual();
        }
    }
    if (correctorName == outerCorrName)
    {
        flowSolutionFinished_ =
            (
                UResidual_.solution < UResidual_.input
             && pResidual_.solution < pResidual_.input
            )
         || (corrector + 1) >= maxOuterCorrectors_;
    }

    const word sedimentCorrName("corrector:" + field_->name());
    if (flowSolutionFinished_ && correctorName == sedimentCorrName)
    {
        fieldCorrectorFinished_ =
            fieldResidual_.solution < fieldResidual_.input
         || corrector >= maxSedimentCorrectors_;
    }

    if
    (
        correctorName == sedimentCorrName
     && flowSolutionFinished_
     && !fieldCorrectorFinished_
    )
    {
        Info<< nl << "Sediment corrector: " << corrector << nl;
    }

    if (correctorName == outerCorrName)
    {
        return flowSolutionFinished_ ? true : false;
    }
    else if (correctorName == sedimentCorrName && !flowSolutionFinished_)
    {
        return true;
    }
    else if (correctorName == sedimentCorrName && flowSolutionFinished_)
    {
        return fieldCorrectorFinished_ ? true : false;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
