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
    (c) 2015 OpenCFD Ltd.
    (c) 2024-2025 Engys Ltd.

Class
    Foam::combustionModels::diffusionMulticomponent

Group
    grpCombustionModels

Description

    Diffusion based turbulent combustion model for multicomponent species.

    The model calculates the laminar finite rate source terms based on
    the kinetic for each reaction in order to begin the combustion and
    evaluates the minimum between this and the cross diffusion rate term
    defined as D*prob*muEff*mag(grad(Yi)*grad(Yj)) if laminarIgn is true.

    where:

    D :     is a model dynamic constant defined as C*f02 where:
            C is a model constant
            f02 = 1  + sqr(O2/oxidantRes), oxidantRes is an user input

    muEff : is the effective turbulent viscosity
    prob :  is a normalized Gaussian shaped distribution around the stoichiometric
            value of each reaction. The distribtion is controled by 'sigma'
            for standard deviation and ftCorr for correction of the stoichiometric
            value.

    In the combustion properties dictionary:

    diffusionMulticomponentCoeffs
    {
        Ci          (1.0 1.0);      // Default to 1
        fuels       (CH4 CO);
        oxidants    (O2 O2);
        YoxStream   (0.23 0.23);    // Default to 0.23
        YfStream    (1.0 1.0);      // Default to 1.0
        sigma       (0.02 0.02);    // Default to 0.02
        oxidantRes  (0.025 0.005);
        ftCorr      (0.0 0.0);      // Default to 0.0
        laminarIgn  false;          // Default false
    }


SourceFiles
    diffusionMulticomponent.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionMulticomponent_H
#define diffusionMulticomponent_H

#include "primitives/Scalar/lists/scalarList.H"
#include "memory/tmp/tmp.H"
#include "reaction/reaction/reaction.H"
#include "memory/autoPtr/autoPtr.H"
#include "combustionModel/combustionModel.H"
#include "basicChemistryModel/basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    Class diffusionMulticomponent Declaration
\*---------------------------------------------------------------------------*/

class diffusionMulticomponent
:
    public combustionModel
{
    // Private data

        //- Run chemistry correction on every outer iteration. Default false.
        bool outerCorrect_;

        //- Reference to the mixture
        const basicSpecieMixture& mixture_;

        //- Reactions
        const PtrList<reaction> reactions_;

        //- Pointer list of source terms
        PtrList<volScalarField> RijPtr_;

        //- Model constants
        scalarList Ci_;

        //- List of fuels for each reaction
        wordList fuelNames_;

        //- List of oxidants for each reaction
        wordList oxidantNames_;

        //- Heat of combustion [J/Kg]
        scalarList qFuel_;

        //- Stoichiometric air-fuel mass ratio
        scalarList stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        scalarList s_;

        //- Oxydaser sream mass concentrations
        scalarList YoxStream_;

        //- Fuel stream mass concentrations
        scalarList YfStream_;

        //- Mean distribution for gaussian probabililty
        scalarList sigma_;

        //- Residual oxydaser
        scalarList oxidantRes_;

        //- ft stochiometric correction
        scalarList ftCorr_;

        //- Relaxatnio factor on total source
        scalar alpha_;

        //- Switch on to laminar combustion for ignition
        bool laminarIgn_;

        //- The time index of the last correction
        label timeIndex_;

        //- Pointer to chemistry model
        autoPtr<basicChemistryModel> chemistryPtr_;


    // Private Member Functions

        //- Initialize
        void init();

        //- Load reaction list
        PtrList<reaction> loadReactions
        (
            const speciesTable& species,
            const dictionary& dict
        ) const;


public:

    //- Runtime type information
    TypeName("diffusionMulticomponent");


    // Constructors

        //- Construct from components
        diffusionMulticomponent
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        diffusionMulticomponent(const diffusionMulticomponent&) = delete;


    //- Destructor
    virtual ~diffusionMulticomponent();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate calculated from fuel consumption rate matrix
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const diffusionMulticomponent&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
