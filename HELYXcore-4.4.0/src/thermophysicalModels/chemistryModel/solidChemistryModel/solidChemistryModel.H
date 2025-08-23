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
    (c) 2013-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::solidChemistryModel

Description
    Extends base solid chemistry model by adding a thermo package, and ODE
    functions.

    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    solidChemistryModelI.H
    solidChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidChemistryModel_H
#define solidChemistryModel_H

#include "reaction/Reactions/ReactionList/ReactionList.H"
#include "ODESystem/ODESystem.H"
#include "fields/volFields/volFields.H"
#include "matrices/simpleMatrix/simpleMatrix.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     Class solidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class SolidThermo>
class solidChemistryModel
:
    public basicSolidChemistryModel,
    public ODESystem
{
protected:

        //- Reference to solid mass fractions
        PtrList<volScalarField>& Ys_;

        //- Thermodynamic data of solids
        const PtrList<SolidThermo>& solidThermo_;

        //- Reactions
        const ReactionList<SolidThermo> reactions_;

        //- Number of solid components
        label nSolids_;

        //- Number of solid reactions
        label nReaction_;

        //- List of reaction rate per solid [kg/m3/s]
        PtrList<volScalarField::Internal> RRs_;

        //- List of active reacting cells
        List<bool> reactingCells_;

        //- Pointer to the Hc models
        matScalarUTable hfModels_;

        //- Pointer to the hs models
        matScalarUTable hsModels_;

        //- Pointer to the rho models
        matScalarUTable rhoModels_;

        //- Pointer to the rho models
        matScalarUTable CpModels_;


    // Protected Member Functions

        //- Write access to source terms for solids
        inline PtrList<volScalarField::Internal>& RRs();

        //- Set reacting status of cell, celli
        void setCellReacting(const label celli, const bool active);


public:

    //- Runtime type information
    TypeName("solidChemistryModel");


    // Constructors

        //- Construct from thermo
        solidChemistryModel(const solidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        solidChemistryModel(const solidChemistryModel&) = delete;


    //- Destructor
    virtual ~solidChemistryModel();


    // Member Functions

        //- The reactions
        inline const ReactionList<SolidThermo>& reactions() const;

        //- The number of reactions
        inline label nReaction() const;


        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            const bool updateC0 = false
        ) const = 0;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<SolidThermo>& r,
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalar& omegaf,
            scalar& omegar
        ) const = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;


        // Solid Chemistry model functions

            //- Return const access to the chemical source terms for solids
            inline const volScalarField::Internal& RRs(const label i) const;

            //- Return total solid source term
            inline tmp<volScalarField::Internal> RRs() const;

            //- Return net solid sensible enthalpy
            inline tmp<DimensionedField<scalar, volMesh>> RRsHs() const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) = 0;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT)
            {
                NotImplemented;
            }

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const
            {
                NotImplemented;
            }

            //- Return the heat release rate [kg/m/s3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const = 0;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const = 0;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const = 0;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidChemistryModel/solidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidChemistryModel/solidChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
