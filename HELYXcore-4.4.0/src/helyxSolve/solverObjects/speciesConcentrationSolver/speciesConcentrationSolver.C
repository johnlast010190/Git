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
    (c) 2019-2025 Engys Ltd.
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "speciesConcentrationSolver.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "multiphaseThermo/multiphaseThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(speciesConcentrationSolver, 0);
}
}

makeFvSolverOption(speciesConcentrationSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::speciesConcentrationSolver::speciesConcentrationSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    turbulencePtr_(nullptr),
    compositionPtr_(nullptr)
{
    Info<< "Reading thermophysical properties\n" << endl;
    //- Phasic and top-level thermo - will point to the same thing if phaseName_
    //  is null
    thermoPtr_ =
        &refCast<rhoThermo>
        (
            multiphaseThermo::lookupOrCreate(obr_, phaseName_)
        );
    globalThermoPtr_ =
        &refCast<rhoThermo>
        (
            basicThermo::lookupOrCreate(obr_, phaseName_)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::speciesConcentrationSolver::read(const dictionary& dict)
{
    // Read the Diffusion coefficients from thermo dict

    if (thermoPtr_ && compositionPtr_)
    {
        Di_.clear();
        Di_.resize(compositionPtr_->Y().size());
        if (thermoPtr_->phaseDict().isDict("speciesDiffusivity"))
        {
            const dictionary& diffDict =
                thermoPtr_->phaseDict().subDict("speciesDiffusivity");
            forAll(Di_, i)
            {
                if (diffDict.found(compositionPtr_->Y()[i].name()))
                {
                    Di_.set
                    (
                        i,
                        Function1<scalar>::New
                        (
                            compositionPtr_->Y()[i].name(),
                            diffDict
                        )
                    );
                }
            }
        }
    }
}


bool Foam::fv::speciesConcentrationSolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    if
    (
        obr_.foundObject<compressible::turbulenceModel>
        (
            compressible::turbulenceModel::propertiesName
        )
    )
    {
        turbulencePtr_ = &obr_.lookupObject<compressible::turbulenceModel>
        (
            compressible::turbulenceModel::propertiesName
        );
    }

    compositionPtr_ = &refCast<fluidMulticomponentThermo>(*thermoPtr_).composition();

    for (auto& Yi : compositionPtr_->Y())
    {
        addField(Yi);

        // For the energy flux added in addSu
        mesh_.schemes().setFluxRequired(Yi.name());
    }

    read(dictionary(dict()));

    return true;
}


void Foam::fv::speciesConcentrationSolver::addField
(
    const volScalarField& field
)
{
    multivariateConvectionFields_.add(field);
}


void Foam::fv::speciesConcentrationSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    const PtrList<volScalarField>& Y = compositionPtr_->Y();
    for (auto& Yi : Y)
    {
        solveNames.append(Yi.name());
        wordList deps;
        deps.append("fvMesh");
        word phiName = "phi";
        if (obr_.foundObject<surfaceScalarField>(phiName) && !isReactingFoam())
        {
            deps.append("U");
        }
        optionalDependencies.insert(Yi.name(), deps);

        correctorMembers.insert
        (
            "nonOrthogonalCorrector:" + Yi.member(),
            {Yi.name()}
        );
    }

    // Make the inert species dependent on all the others as we need them
    // to solve first
    if (compositionPtr_->defaultSpecie() != -1)
    {
        DynamicList<word> activeSpecies;
        forAll(Y, speciesI)
        {
            if (speciesI != compositionPtr_->defaultSpecie())
            {
                activeSpecies.append(Y[speciesI].name());
            }
        }
        requiredDependencies.set
        (
            Y[compositionPtr_->defaultSpecie()].name(),
            activeSpecies
        );
    }

    optionalDependencies.insert("turbulence", solveNames);

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


bool Foam::fv::speciesConcentrationSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        return (corrector >= solnControlPtr_->nOuterCorr()-1);
    }
    else if
    (
        correctorName.substr(0, correctorName.find(':'))
     == "nonOrthogonalCorrector"
    )
    {
        // Non-orthogonal correctors
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else
    {
        return false;
    }
}


Foam::fv::convectionScheme<Foam::scalar>&
Foam::fv::speciesConcentrationSolver::convTerm()
{
    const word phiName("phi");
    if
    (
        convectionScheme_.empty()
     && obr_.foundObject<surfaceScalarField>(phiName)
    )
    {
        const word YiName(addPhaseName("Yi"));
        const word schemeName("div(" + phiName + "," + YiName + ")");
        convectionScheme_ =
            fv::convectionScheme<scalar>::New
            (
                mesh_,
                multivariateConvectionFields_,
                obr_.lookupObjectRef<surfaceScalarField>(phiName),
                mesh_.schemes().divScheme(schemeName)
            );
    }
    return convectionScheme_.ref();
}


tmp<fvScalarMatrix>
Foam::fv::speciesConcentrationSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    volScalarField& rho =
        obr_.lookupObjectRef<volScalarField>(addPhaseName("rho"));

    // If phi not found (no flow solver running), we include no convection
    word phiName = "phi";
    word YiName = addPhaseName("Yi");

    PtrList<volScalarField>& Y = compositionPtr_->Y();
    forAll(Y, specieI)
    {
        if (Y[specieI].name() == fieldName)
        {
            if (compositionPtr_->active(specieI))
            {
                volScalarField& Yi = Y[specieI];

                fv::options& fvOptions = this->fvOptions();

                // Store diffusion coefficient for use by coupling BCs
                const word rhoDName = "rhoD_" + Yi.name();
                if (!obr_.foundObject<volScalarField>(rhoDName))
                {
                    tmp<volScalarField> rhoD
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                rhoDName,
                                mesh_.time().timeName(),
                                obr_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh_,
                            dimDensity*sqr(dimLength)/dimTime
                        )
                    );
                    regIOobject::store(rhoD.ptr());
                }
                volScalarField& rhoD =
                    obr_.lookupObjectRef<volScalarField>(rhoDName);

                // Laminar diffusion coefficient
                // Not (yet) integrated into thermo library as this depends on
                // the diffusion model (Fick vs multi-component, etc)
                if (Di_.set(specieI))
                {
                    // Specified diffusion vs temp
                    volScalarField D
                    (
                        IOobject
                        (
                            "D_" + YiName,
                            mesh_.time().timeName(),
                            obr_
                        ),
                        mesh_,
                        sqr(dimLength)/dimTime,
                        zeroGradientFvPatchScalarField::typeName
                    );
                    D.primitiveFieldRef() =
                        Di_[specieI].value(thermoPtr_->T().primitiveField())();
                    D.correctBoundaryConditions();

                    rhoD = rho*D;
                }
                else
                {
                    // If not specified, default to laminar viscosity (Reynolds
                    // analogy)
                    rhoD = thermoPtr_->mu();
                }

                // Turbulent contribution
                if (turbulencePtr_)
                {
                    rhoD += turbulencePtr_->Dmt();
                }

                // Tortuosity adjustment
                if (thermoPtr_->properties().found("tortuosity"))
                {
                    scalar tortuosity =
                        thermoPtr_->properties().lookup<scalar>("tortuosity");
                    rhoD /= tortuosity;
                }

                if (compositionPtr_->solve(specieI))
                {
                    tmp<fvScalarMatrix> tYiEqn
                    (
                        fvm::ddt(rho, Yi, "ddt(" + rho.name() + ",Yi)")
                      - fvm::laplacian
                        (
                            rhoD, Yi, "laplacian(rhoD_Yi,Yi)"
                        )
                    ==
                        fvOptions(rho, Yi)
                    );

                    const surfaceScalarField* phiPtr =
                        obr_.lookupObjectRefPtr<surfaceScalarField>(phiName);
                    if (phiPtr)
                    {
                        tYiEqn.ref() += convTerm().fvmDiv(*phiPtr, Yi);
                    }

                    if (thermoPtr_->properties().found("porosity"))
                    {
                        scalar eps =
                            thermoPtr_->properties().lookup<scalar>("porosity");
                        tYiEqn.ref() *= eps;
                    }

                    tYiEqn.ref().relax(YiName);
                    fvOptions.constrain(tYiEqn.ref());

                    dictName = YiName;
                    return tYiEqn;
                }
            }
        }
    }
    return tmp<fvScalarMatrix>();
}


void Foam::fv::speciesConcentrationSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    forAll(Y, specieI)
    {
        if (Y[specieI].name() == solveName)
        {
            if (compositionPtr_->solve(specieI))
            {
                this->fvOptions().correct(Y[specieI]);
            }
            else if (specieI == compositionPtr_->defaultSpecie())
            {
                // This will switch the specie to calculated type of bc
                // and also ensures that all boundary values are 0
                compositionPtr_->normalise();

                if (!isReactingFoam())
                {
                    // Inert species is the last to be solved so we
                    // wait until here to apply the density correction
                    globalThermoPtr_->correct();

                    // For stability, keep the solver density updated
                    volScalarField& rho =
                        obr_.lookupObjectRef<volScalarField>("rho");
                    rho = globalThermoPtr_->rho();
                }
            }
            break;
        }
    }
}


void Foam::fv::speciesConcentrationSolver::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (correctorName == solverObject::outerCorrectorName)
    {
        if (convectionScheme_.valid())
        {
            convectionScheme_.clear();
        }
    }
}


void Foam::fv::speciesConcentrationSolver::getSourceGraph
(
    wordList& fieldNames,
    HashTable<wordList>& sourceDependencies
)
{
    // We act as a source for the energy equation, to add the energy flux
    // due to mass transfer

    // Equations for which we are providing a source term
    fieldNames = {globalThermoPtr_->he().name()};

    // Require the mass fractions to have been solved before these terms are
    // added to energy equation
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    DynamicList<word> speciesNames;
    forAll(Y, specieI)
    {
        speciesNames.append(Y[specieI].name());
    }
    sourceDependencies.insert(fieldNames[0], speciesNames);
}


void Foam::fv::speciesConcentrationSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (isReactingFoam())
    {
        return;
    }

    // Adding terms to the energy equation to account for energy flux
    // due to mass transfer
    PtrList<volScalarField>& Y = compositionPtr_->Y();
    const volScalarField& he = globalThermoPtr_->he();
    volScalarField hej("hej", he);

    forAll(Y, specieI)
    {
        if (compositionPtr_->active(specieI))
        {
            volScalarField& Yi = Y[specieI];

            const word rhoDName = "rhoD_" + Yi.name();
            volScalarField& rhoD =
                obr_.lookupObjectRef<volScalarField>(rhoDName);

            const objectRegistry& matObr = obr_.subRegistry("materialModels");
            const materialTables& matTable =
                matObr.lookupObject<materialTables>("materialTables");
            hej.forceAssign
            (
                matTable(heModel::typeName, phaseName_, Yi.name())()
            );

            // Species diffusive mass flux
            tmp<surfaceScalarField> tJi =
               -fvm::laplacian(rhoD, Yi, "laplacian(rhoD_Yi,Yi)")->flux();
            const surfaceScalarField& Ji = tJi();

            tmp<fv::convectionScheme<scalar>> convScheme =
                fv::convectionScheme<scalar>::New
                (
                    mesh_,
                    Ji,
                    mesh_.schemes().divScheme
                    (
                        "div(Ji,"
                      + IOobject::groupName
                        (
                            eqn.psi().member() + "i",
                            eqn.psi().group()
                        )
                      + ")"
                    )
                );

            // Clip it to small value to avoid devision by zero.
            // (this could be problematic for single precision)
            tmp<surfaceScalarField> interpolateJiHe
            (
                convScheme->interpolate(Ji, eqn.psi())
            );
            forAll(interpolateJiHe.ref(), facei)
            {
                scalar& value = interpolateJiHe.ref()[facei];
                if (mag(value) < SMALL)
                {
                    value = SMALL;
                }
            }
            forAll(interpolateJiHe.ref().boundaryField(), patchi)
            {
                fvsPatchField<scalar>& bc =
                    interpolateJiHe.ref().boundaryFieldRef()[patchi];
                if (!bc.fixesValue())
                {
                    scalarField& values = dynamic_cast<scalarField&>(bc);
                    forAll(values, patchFacei)
                    {
                        if (mag(values[patchFacei]) < SMALL)
                        {
                            values[patchFacei] = SMALL;
                        }
                    }
                }
            }

            // We want convection of the species energy, div(Ji, hej), but
            // treating implicitly through he.
            // eqn is RHS, so subtract
            eqn -=
                convScheme->fvmDiv
                (
                    Ji*convScheme->interpolate(Ji, hej)/interpolateJiHe,
                    eqn.psi()
                );
        }
    }
}


// ************************************************************************* //
