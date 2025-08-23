/*---------------------------------------------------------------------------* \
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

Class
    Foam::fv::multiphaseSegregatedFlowSolver

Group
    grpFvOptionsSolvers

Description

    solverOption for fluid pressure-based U-p solve

SourceFiles
    multiphaseSegregatedFlowSolver.C

\*---------------------------------------------------------------------------*/
#ifndef multiphaseSegregatedFlowSolver_H
#define multiphaseSegregatedFlowSolver_H

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "rhoThermo/rhoThermo.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "cfdTools/general/pressureControl/pressureControl.H"
#include "solverObjects/flowSolver/flowSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiphaseThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class multiphaseSegregatedFlowSolver Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSegregatedFlowSolver
:
    public flowSolver
{
protected:
    volScalarField* p_;
    autoPtr<volVectorField> U_;
    multiphaseThermo *thermoPtr_;
    autoPtr<volScalarField> rho_;

    //- Volumetric flux, relative to mesh motion and MRF
    autoPtr<surfaceScalarField> phiv_;

    //- Mass flux, relative to mesh motion and MRF
    autoPtr<surfaceScalarField> phi_;

    autoPtr<surfaceVectorField> Uf_;
    bool closedVolume_;
    autoPtr<dimensionedScalar> initialMass_;
    autoPtr<dimensionedScalar> compressibility_;
    autoPtr<dimensionedScalar> rhoMin_;
    autoPtr<dimensionedScalar> rhoMax_;
    compressible::turbulenceModel* turbulence_;
    Switch printContErr_;
    scalar cumulativeContErr_;
    scalar maxCo_;

    autoPtr<Function1<scalar>> maxCoDataPtr_;
    solutionControl* solnControlPtr_;
    autoPtr<pressureControl> pressureControl_;

    //- Is the case buoyant?
    bool buoyant_;

    //- Density used for buoyancy term
    autoPtr<volScalarField> bRho_;

    //- Indicates a correctPhi is pending, due to topo change or initial init
    bool needsCorrectPhi_;

    //- Used to flag initialisation of phiv and phi which needs to be delayed
    //  until the first mesh move
    bool phiNeedsInit_;

    //- Initial velocity solve is implicit
    bool implicitUSolve_;
    //- Velocity corrector solves are implicit
    bool implicitUCorrector_;
    //- Use a momentum predictor
    bool momentumPredictor_;

    label corr_;
    label nonOrthCorr_;

    tmp<volVectorField> HbyA_;
    tmp<surfaceScalarField> phivHbyA_;
    tmp<volScalarField> rAU_;
    tmp<surfaceScalarField> rAUf_;
    tmp<volScalarField> p0_;
    tmp<volScalarField> divU_;

    // Use tmps for proper reference counting (allows cleanup from outside)
    tmp<fvVectorMatrix> UEqnLHS_;
    tmp<fvScalarMatrix> pDDtEqn_;
    tmp<fvScalarMatrix> pEqn_;

    tmp<fvScalarMatrix> pEqnIncomp_;
    PtrList<fvScalarMatrix> pEqnComps_;

    //- Name of the rAUf field
    static const word rAUfName_;

protected:
    void assembleUEqnLHS();
    void assemblepEqn();
    void solveContinuity();
    void printIncoContinuityErrors(const word&);
    void printCompContinuityErrors(const word&);

    //- Return pressure
    const volScalarField& p() const override
    {
        return *p_;
    }

    //- Return thermo model
    const fluidThermo& thermo() const override
    {
        NotImplemented;
    }

    //- Return solver rho
    tmp<volScalarField> rho() const override
    {
        return tmp<volScalarField>(rho_());
    }

    tmp<volScalarField> buoyantRho() const override
    {
        return rho();
    }

    const solutionControl& solnControl() const override
    {
        return *solnControlPtr_;
    }

    const surfaceScalarField& phiv() const override
    {
        return phiv_();
    }

    void correctPhi(bool defaultCorrectPhi);

public:

//- Runtime type information
    TypeName("multiphaseSegregatedFlowSolver");

// Constructors

    // - Construct from components
    multiphaseSegregatedFlowSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );


//- Destructor
    virtual ~multiphaseSegregatedFlowSolver()
    {}

// Member Functions

    virtual bool initialise() override;

    virtual void getSolveGraph
    (
        wordList& solveNames,
        HashTable<wordList>& requiredDependencies,
        HashTable<wordList>& optionalDependencies,
        HashTable<wordList>& correctorMembers
    ) override;

    virtual scalar getMaxTimeStep() override;

    virtual bool isFinalCorrector
    (
        const label corrector,
        const word& correctorName
    ) override;

    virtual void beginIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

    virtual tmp<surfaceScalarField>
    faceBuoyancyForce
    (
        const word& pDiffName, bool includeSnGradP=false
    ) const override;

    virtual tmp<fvScalarMatrix> assembleScalarMatrix
    (
        const word& fieldName,
        bool& finalSolve,
        word& dictName
    ) override;
    virtual tmp<fvVectorMatrix> assembleVectorMatrix
    (
        const word& fieldName,
        bool& finalSolve,
        word& dictName
    ) override;

    virtual void correct(const word&, const word&) override;

    virtual void endIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

    virtual bool movePoints() override;

    virtual void topoChange(const polyTopoChangeMap& map) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
