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
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fv::porousSolidThermalEquilibrium

Group
    grpFvOptionsSolvers

Description
    Calculates the solid energy assuming thermal equilibrium in a porous
    matrix. Adds required source terms to fluid energy equation to simulate
    effective properties. Only formally correct where porosity is uniform.

SourceFiles
    porousSolidThermalEquilibrium.C

SeeAlso
    Foam::fv::fluidEnergySolver

\*---------------------------------------------------------------------------*/

#ifndef porousSolidThermalEquilibrium_H
#define porousSolidThermalEquilibrium_H

#include "solverObjects/solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "basicThermo/basicThermo.H"
#include "combustionModel/combustionModel.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvSolutionRegistry;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class porousSolidThermalEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class porousSolidThermalEquilibrium
:
    public solverObject
{
protected:

    //- The fluid reaction
    combustionModel* reactionPtr_;

    //- The thermos
    basicThermo* thermoPtr_;
    basicThermo* fluidThermoPtr_;

    word fluidRegion_;
    const fvSolutionRegistry& fluidDb_;

    //- Correction to the effective density - stored for time derivative
    autoPtr<volScalarField> porousRhoEffCorr_;

    solutionControl* solnControlPtr_;

    //- State variable to avoid recursively collecting source dependencies from
    //  outself
    bool addingDependencies_;


public:

    //- Runtime type information
    TypeName("porousSolidThermalEquilibrium");

    // Constructors

        // - Construct from components
        porousSolidThermalEquilibrium
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Dissallow default bitwise copy construction
        porousSolidThermalEquilibrium(porousSolidThermalEquilibrium&) = delete;


    //- Destructor
    virtual ~porousSolidThermalEquilibrium()
    {}


    // Member Functions

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& derivedSolves,
            HashTable<solveList>& requiredDependencies,
            HashTable<solveList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;

        virtual void getSourceGraph
        (
            solveList& fields,
            SolveTable<solveList>& sourceDependencies
        ) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

        //- Read the porousSolidThermalEquilibrium data
        virtual void read(const dictionary& dict) override;


    // Member Operators

        //- Dissallow default bitwise assignment
        void operator=(porousSolidThermalEquilibrium&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
