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
    (c) 2025 Engys Ltd.

Class
    Foam::fv::pcweAcousticSolver

Group
    grpFvOptionsSolvers

Description
    Solves the Perturbed Convective Wave Equation (PCWE) for acoustics:

    \f[
    \frac{1}{c_0^2} \frac{D^2 \psi_a}{Dt^2} - \nabla^2 \psi_a = -\frac{1}{\rho_0 c_0^2} \frac{D p_{ic}}{Dt}
    \f]

    where:
    - \f$ \psi_a \f$ is the acoustic potential
    - \f$ p_{ic} \f$ is the incompressible flow pressure
    - \f$ c_0 \f$ is the speed of sound

Usage
    \verbatim

        acousticSolver
        {
            type             pcweAcousticSolver;
            convectiveTerms  no;
            flowSolverSource yes;
            startTime        0;
            dtAco            1e-06;
        }

    \endverbatim

    Where:

    \table
        Property          | Description                               | Required | Default value
        type              | Type of the acoustic solver               | yes      | pcweAcousticSolver
        convectiveTerms   | Include convective terms (yes/no)         | no       | no
        flowSolverSource  | Use source from the flow solver (yes/no)  | no       | no
        startTime         | Start time of the acoustic solver         | no       | 0
        dtAco             | Acoustic time step                        | no       | time-step value controlDict
    \endtable

SourceFiles
    pcweAcousticSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pcweAcousticSolver_H
#define pcweAcousticSolver_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class pcweAcousticSolver Declaration
\*---------------------------------------------------------------------------*/

class pcweAcousticSolver
:
    public solverObject
{
protected:

    solutionControl* solnControlPtr_;

    basicThermo* thermoPtr_;
    surfaceScalarField* phiPtr_;
    const volVectorField* UPtr_;
    const volScalarField* pPtr_;

    //- Speed of sound [m/s]
    tmp<volScalarField> c0_;

    //- Start time to solve the PCWE [s]
    scalar startTime_;

    //- Acoustic time step [s]
    scalar dtAco_;

    //- Acoustic potential field
    volScalarField psiAco_;

    //- Acoustic pressure field
    volScalarField pAco_;

    //- Total derivative of the incompressible pressure field
    //  equal to time derivative if convection is not included
    volScalarField dpdt_;

    //- Flag to indicate that the convective terms are included in the PCWE
    bool convectiveTerms_;

    //- Flag to indicate if the source term is from the pressure of flow solver
    bool flowSolverSource_;

    //- Number of sub-cycles to solve the PCWE
    label nAcoSubCycles_;

public:

    //- Runtime type information
    TypeName("pcweAcousticSolver");

    // Constructors

        // - Construct from components
        pcweAcousticSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        pcweAcousticSolver(pcweAcousticSolver&) = delete;

        //- Destructor
        virtual ~pcweAcousticSolver()
        {}


    // Member Functions

        virtual void read(const dictionary& dict) override;

        //- Initialize solver (fields/variables)
        virtual bool initialise() override;

        //- Identify any solves handled by this solverObject,
        // their levels of iteration and dependent solves
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

        //- Assemble and solve inisde separate function since
        //  optional sub-cycling otherwise duplicates all code
        virtual void solveEquation
        (
            volScalarField& field
        );

        virtual void correct(const word&, const word&) override;

        //- Compute the number of sub-cycles needed to satisfy the
        //  target time step for acoustics
        virtual label nSubCycles() const;

    // Member Operators

        //- Dissallow default bitwise copy assignment
        void operator=(pcweAcousticSolver&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
