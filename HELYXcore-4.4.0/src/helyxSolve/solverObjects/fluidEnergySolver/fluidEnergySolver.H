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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::fluidEnergySolver

Group
    grpFvOptionsSolvers

Description

    solverObject for enthalpy or energy solution with optional porous solid
    (thermal equilibrium assumed)

SourceFiles
    fluidEnergySolver.C

\*---------------------------------------------------------------------------*/
#ifndef fluidEnergySolver_H
#define fluidEnergySolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "primitives/enums/Enum.H"
#include "eulerianMultiphaseSystem/eulerianMultiphaseSystem.H"
#include "solverObjects/speciesConcentrationSolver/speciesConcentrationSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class fluidEnergySolver Declaration
\*---------------------------------------------------------------------------*/

class fluidEnergySolver
:
    public solverObject
{
    // Private Data Members

        //- Convenience typedef
        typedef speciesConcentrationSolver speciesType;


protected:

    const solutionControl* solnControlPtr_;

    basicThermo* thermoPtr_;

    //- Pointer to phase system for multiple phasic energy models
    eulerianMultiphaseSystem* phaseSystemPtr_;

    //- Does it have MRF?
    bool hasMRF_;

    //- Stabilisation mechanism for low Froude number flows
    Switch lowFroudeStabilisation_;

    //- Coefficient for stabilisation of low Froude number flows (0-1)
    scalar lowFroudeStabilisationCoeff_;

    //- Switch to add or not the viscous term
    Switch viscousStressTerm_;

    //- Number of energy correctors for phasic system
    label nEnergyCorrectors_;

    enum formulationType
    {
        totalEnergyEnergy = 0,
        totalEnergyTemperature,
        boussinesqEnergy,
        boussinesqTemperature,
        phasicEnergy
    } formulation_;

    autoPtr<volScalarField> K_;

    autoPtr<eulerianPhaseSystem::heatTransferTable> heatTransferPtr_;

    //- Create K field
    virtual void createK();

    //- Calculate MRF flux
    tmp<volScalarField> divMRF
    (
        const surfaceScalarField& phi,
        const volScalarField* alpha = nullptr
    );

    //- Calculate viscous stress term
    tmp<volScalarField> viscousStressTerm(const volVectorField&) const;


public:

    static const Enum<formulationType> formulationTypeNames_;

    //- Runtime type information
    TypeName("fluidEnergySolver");


    // Constructors

        // - Construct from components
        fluidEnergySolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fluidEnergySolver()
    {}


    // Member Functions

        //- Helper function to handle multivariable convection scheme
        tmp<fvScalarMatrix> fvmDiv
        (
            const surfaceScalarField& phi,
            const volScalarField& U
        );

        //- Add energy to multivariable scheme
        void addEnergyToMultivariableScheme();

        //- Add energy to multivariable scheme
        tmp<volScalarField> lookupOrConstructDpdt
        (
            const volScalarField& p
        );

        //- Initialise solver object
        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual bool isFinalCorrector
        (
            const label corrector,
            const word& correctorName
        ) override;

        virtual tmp<fvScalarMatrix> assembleScalarMatrix
        (
            const word& solveName,
            bool& finalSolve,
            word& dictName
        ) override;

        virtual tmp<fvScalarMatrix> assemblePhasicEEqn(const label phasei);

        virtual void correct
        (
            const word& solveName,
            const word& regionName
        ) override;

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
