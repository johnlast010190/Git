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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::regionModels::reactingOneDim

Description
    Reacting, 1-D pyrolysis model
    NOTE: The moveMesh option can only be applied to solid reaction such as:
          PMMA -> gas at the moment.

SourceFiles
    reactingOneDim.C

\*---------------------------------------------------------------------------*/

#ifndef reactingOneDim_H
#define reactingOneDim_H

#include "pyrolysisModels/pyrolysisModel/pyrolysisModel.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "radiationModels/radiationModel/radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class reactingOneDim Declaration
\*---------------------------------------------------------------------------*/

class reactingOneDim
:
    public pyrolysisModel
{
    // Private member functions

        //- Read model controls
        void readReactingOneDimControls();


protected:

        // Reference to solid thermo properties

            //- List of solid components
            PtrList<volScalarField>& Ys_;


        // Solution parameters

            //- Minimum delta for combustion
            scalar minimumDelta_;


        // Fields

            //- Total gas mass flux to the primary region [kg/m2/s]
            surfaceScalarField phiGas_;

            //- Sensible enthalpy gas flux [J/m2/s]
            volScalarField phiHsGas_;

            //- Heat release rate [J/s/m3]
            volScalarField chemistryQdot_;


        // Source term fields

            //- In depth radiative heat flux [W/m2]
            volScalarField qr_;


        // Checks

            //- Cumulative lost mass of the condensed phase [kg]
            dimensionedScalar lostSolidMass_;

            //- Cumulative mass generation of the gas phase [kg]
            dimensionedScalar addedGasMass_;

            //- Total mass gas flux at the pyrolysing walls [kg/s]
            scalar totalGasMassFlux_;

            //- Total heat release rate [J/s]
            dimensionedScalar totalHeatRR_;


        // Options

            //- Add in depth radiation source term
            bool  qrHSource_;

            //- Use chemistry solvers (ode or sequential)
            bool useChemistrySolvers_;


    // Protected member functions

        //- Read control parameters from dictionary
        bool read() override;

        //- Read control parameters from dict
        bool read(const dictionary& dict) override;

        //- Update/move mesh based on change in mass
        void topoChange(const scalarField& mass0);

        //- Update radiative flux in pyrolysis region
        void updateqr();

        //- Mass check
        void calculateMassTransfer();

        //- Solve continuity equation
        void solveContinuity();

        //- Solve energy
        void solveEnergy();

        //- Solve solid species mass conservation
        void solveSpeciesMass();


public:

    //- Runtime type information
    TypeName("reactingOneDim");


    // Constructors

        //- Construct from type name and mesh
        reactingOneDim
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name, mesh and dictionary
        reactingOneDim
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );

        //- Disallow default bitwise copy construction
        reactingOneDim(const reactingOneDim&) = delete;


    //- Destructor
    virtual ~reactingOneDim();


    // Member Functions

        //- Return the total gas mass flux to primary region [kg/m2/s]
        virtual const surfaceScalarField& phiGas() const override;

        //- External hook to add mass to the primary region
        //  (patchi on primary region, facei of patchi)
        virtual scalar addMassSources
        (
            const label patchi,
            const label facei
        ) override;

        //- Mean diffusion number of the solid region
        virtual scalar solidRegionDiffNo() const override;

        //- Evolve the pyrolysis equations
        virtual void evolveRegion() override;

        //- Provide some feedback
        virtual void info() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactingOneDim&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
