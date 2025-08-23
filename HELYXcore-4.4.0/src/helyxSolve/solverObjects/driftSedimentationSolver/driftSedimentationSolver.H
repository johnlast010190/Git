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
    (c) 2011 Chenxiaoxiao
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fv::driftSedimentationSolver

Group
    grpFvOptionsSolvers

Description
    Drift sedimentation solver for sedimentation.

Reference:
    Chen, X., & Yu, Z. (2022). DriftScalarDyFoam: An OpenFOAM-Based Multistage
    Solver for Drifting Snow and Its Distribution Around Buildings. Frontiers
    in Earth Science, 10, 822140. https://doi.org/10.3389/feart.2022.822140

SourceFiles
    driftSedimentationSolver.C

\*---------------------------------------------------------------------------*/
#ifndef driftSedimentationSolver_H
#define driftSedimentationSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "rhoThermo/rhoThermo.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"

#include "materialModels/materialTables/materialTables.H"
#include "materialModels/baseModels/materialModels.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
/*---------------------------------------------------------------------------*\
                    Class driftSedimentationSolver Declaration
\*---------------------------------------------------------------------------*/

class driftSedimentationSolver
:
    public solverObject
{
protected:

    // Protected Data

        //- Pointer to solution controls
        const solutionControl* solnControlPtr_;

        //- Field name
        word fieldName_;

        //- To store field
        autoPtr<volScalarField> field_;

        //- Falling velocity sediment velocity
        autoPtr<volVectorField> wf_;

        //- Flux of falling sediment
        autoPtr<surfaceScalarField> phiWf_;

        //- Total mass exchange
        autoPtr<volScalarField> M_;

        //- Sediment accumulated height
        autoPtr<volScalarField> deltaH_;

        //- Sediment accumulated height
        autoPtr<volScalarField> deltaHold_;

        //- Surface level
        autoPtr<volVectorField> surfaceLevel_;

        //- Depth of the sediment [m]
        autoPtr<volScalarField> totalSedimentDepth_;

        //- Sediment depth per unit area [kg/m^2]
        autoPtr<volScalarField> sedimentLoad_;

        //- Wall shear stress
        autoPtr<volVectorField> shearStress_;

        //- Sediment patches
        labelList sedimentPatchIDs_;

        //- Errosion constant
        scalar erosionConstant_;

        //- Threshold shear wind speed
        scalar thresholdShearWindSpeed_;

        //- List of locations of cells with 0.1 distance from
        //  the sedimentation patches
        HashTable<autoPtr<labelList>, label> patchesCells_;

        //- Distance for U internal contribution to the sedimentation
        scalar UInternalDistance_;

        //- Patches to be included for U internal sedimentation
        boolList UInternalPatches_;

        //- Residual control object
        struct residualControl
        {
            scalar input;
            scalar solution;
        };

        //- Velocity residual
        residualControl UResidual_;

        //- Pressure residual
        residualControl pResidual_;

        //- Field concentration residual
        residualControl fieldResidual_;

        //- Max number of outer correctors
        label maxOuterCorrectors_;

        //- Max number of outer correctors
        label maxSedimentCorrectors_;

        //- Sediment density
        scalar rhoSediment_;

        //- Flow solution finished?
        bool flowSolutionFinished_;

        //- Field at least ones solved?
        bool fieldSolved_;

        //- If phasic - the phasic thermo
        rhoThermo* thermoPtr_;

        //- Pointer to turbulence model
        const compressible::turbulenceModel* turbulencePtr_;

        //- Normal vector (opossite to gravity direction)
        vector normal_;

        //- Starting surface name
        word startingSurfaceName_;

        //- Is the filed solution finished?
        bool fieldCorrectorFinished_;

        //- Alternative errosion equation
        Switch alternativeErosion_;

        //- Tominaga deposition rate
        Switch TominagaDepositionRate_;

        //- Is deposition added when erosion ocours?
        Switch depositionWithErosion_;


public:

    //- Runtime type information
    TypeName("driftSedimentationSolver");


    // Constructors

        // - Construct from components
        driftSedimentationSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

    //- Destructor
    virtual ~driftSedimentationSolver()
    {}


    // Member Functions

        //- Update phiWf
        void updatePhiWf(bool updateOnlyBoundary = false);

        //- Initialise the object (useful since the fields are usually already
        //  constructed)
        virtual bool initialise() override;

        //- Identify any solves handled by this solverObject. These can be
        //  any unique identifier of each of the solution processes;
        //  optionally the name of a field, but not necessarily. For each
        //  solve, optionally also record the names of any additional fields
        //  that are derived from the solve. Also return:
        //   * Any dependent solves (solve names or fields) which must be
        //     completed first.
        //   * The correctors which this solver object is part of.
        //   * For each corrector, return the solves which should be part of
        //     it. If a corrector name is given, all solves that are part
        //     of that corrector will be added.
        //  The supplied tables are guaranteed empty when called.
        //  Only override one of the three getSolveGraph() signatures below.
        //  The first allows for dependent fields in different regions to be
        //  specified, as well as derived fields; the second allows for
        //  derived fields only; the third allows neither.
        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        //- Called at the top of the corrector loop to determine if
        //  enough correctors have been run at this level (including
        //  the one about to be run, which is always run - possibly the
        //  final corrector). Correctors will be continued until
        //  all solvers have returned true here
        virtual bool isFinalCorrector
        (
            const label corrector,
            const word& correctorName
        ) override;

        //- Called to return a matrix. Can return nothing, in which
        // case nothing will be solved.
        // finalSolve is initialised with the finalIteration data from the
        // mesh. Overwrite it to override the 'final solve' selection passed
        // to the solver. dictName is initialised with the word::null and
        // can be overwritten to specify the solver name to be selected in
        // fvSolution.
        virtual tmp<fvScalarMatrix> assembleScalarMatrix
        (
            const word& fieldName,
            bool& finalSolve,
            word& solveName
        ) override;

        //- Corrector step after solve
        virtual void correct(const word&, const word&) override;

        //- Called at the bottom of the corrector loop
        virtual void endIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        ) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
