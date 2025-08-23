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
    (c) 2023-2025 Engys Ltd.

Class
    Foam::fv::MULESVolumeFractionSolver

Group
    grpFvOptionsSolvers

Description
    Solver object using MULES to transport the volume fraction (two- or multi-
    phase)

SourceFiles
    MULESVolumeFractionSolver.C

\*---------------------------------------------------------------------------*/
#ifndef MULESVolumeFractionSolver_H
#define MULESVolumeFractionSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "coordinate/systems/coordinateSystem.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "materialModels/materialTables/materialTables.H"
#include "materialModels/baseModels/materialModels.H"
#include "eulerianPhaseSystems/eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiphaseThermo;
class solutionControl;
class eulerianMultiphaseSystem;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class MULESVolumeFractionSolver Declaration
\*---------------------------------------------------------------------------*/

class MULESVolumeFractionSolver
:
    public solverObject
{
protected:

    //- Mulitiphase termo model
    multiphaseThermo* thermoPtr_;

    //- Multiphase system (present for Euler-Euler)
    eulerianMultiphaseSystem* phaseSystemPtr_;

    //- Number of alpha-equation sub-cycles
    label nAlphaSubCycles_;

    //- Number of alpha correctors
    label nAlphaCorr_;

    //- Implicit MULES
    bool MULESCorr_;

    //- Whether to apply previous iteration's corrector
    //  Apply the compression correction from the previous iteration
    //  Improves efficiency for steady-simulations but can only be applied
    //  once the alpha field is reasonably steady, i.e. fully developed
    bool alphaApplyPrevCorr_;

    //- Index of the phase which isn't solved for
    label passiveIndex_;

    //- Pointer to volume flux
    const surfaceScalarField* phivPtr_;

    //- Pointer to mass flux (can be null if not present)
    surfaceScalarField* phiPtr_;

    //- Pointer to velocity field (can be null in Euler-Euler case)
    const volVectorField* UPtr_;

    //- Maximum interface Courant number
    scalar maxAlphaCo_;
    autoPtr<Function1<scalar>> maxAlphaCoDataPtr_;

    solutionControl* solnControlPtr_;

    //- Whether alpha scheme is Crank-Nicolson
    bool crankNicolson_;

    //- Crank-Nicolson off-centering coefficient
    scalar cnCoeff_;

    //- Phase fluxes. Only persistent if using Crank Nicolson or phi not present
    PtrList<surfaceScalarField> alphaPhivs_;

    //- Previous iteration correctors
    PtrList<surfaceScalarField> alphaPhivCorr0_;

    //- Local timestepping
    bool LTS_;

    //- Optional phase-pair specific compression coefficients
    autoPtr<HashTable<scalar>> cAlphas_;

    //- Global interface compression
    scalar cAlphaGlobal_;

    //- Optional flag to disable outlet flux correction
    bool correctInflowOutflow_;

    void updateTimeSchemeInfo();

    void calcSuSp
    (
        const label phasei,
        autoPtr<volScalarField::Internal>& Su,
        autoPtr<volScalarField::Internal>& Sp
    );

    virtual void solveAlphas(const volScalarField::Internal& divU);

    void correctInflowOutflow
    (
        surfaceScalarField& alphaPhi,
        const surfaceScalarField& phiv,
        const volScalarField& alpha
    ) const;

    bool isDeprecatedAlphaControl
    (
        const dictionary& alphaControls, const word& key
    );


public:

    //- Runtime type information
    TypeName("MULESVolumeFractionSolver");


    // Constructors

        // - Construct from components
        MULESVolumeFractionSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

    //- Destructor
    virtual ~MULESVolumeFractionSolver()
    {}


    // Member Functions

        virtual void read(const dictionary& dict) override;

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& derivedFields,
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

        virtual void correct(const word&, const word&) override;

        virtual void endIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        ) override;

        virtual void topoChange(const polyTopoChangeMap& map) override;

        //- Set local timestep limits for this function object
        virtual bool setRDeltaT() override;

        //- Is it last time step corrector (sweep/smooth)
        virtual bool isLastRDeltaTCorrector() const override
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
