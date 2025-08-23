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
    (c) 2016-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::RASModels::kOmegaSSTLM

Group
    grpLESTurbulence

Description
    Langtry-Menter 4-equation transitional SST model
    based on the k-omega-SST RAS model. It includes the
    crossflow correction model (optional)

    References:
    \verbatim
        Langtry, R. (2015). Extending the gamma-Retheta
        correlation based transition model for crossflow effects.
        In 45th AIAA fluid dynamics conference.

        Langtry, R. B., & Menter, F. R. (2009).
        Correlation-based transition modeling for unstructured parallelized
        computational fluid dynamics codes.
        AIAA journal, 47(12), 2894-2906.

        Menter, F. R., Langtry, R., & Volker, S. (2006).
        Transition modelling for general purpose CFD codes.
        Flow, turbulence and combustion, 77(1-4), 277-303.

        Langtry, R. B. (2006).
        A correlation-based transition model using local variables for
        unstructured parallelized CFD codes.
        Phd. Thesis, Universität Stuttgart.
    \endverbatim

    The model coefficients are
    \verbatim
        kOmegaSSTCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1;
            c1          10;
            F3          no;

            // Default LM coefficients
            ca1         2;
            ca2         0.06;
            ce1         1;
            ce2         50;
            cThetat     0.03;
            sigmaThetat 2;

            lambdaErr   1e-6;
            maxLambdaIter 30;

            // Crossflow correction
            crossflow true;
            roughnessHeight 3.3e-6; //default 3e-6
        }
    \endverbatim

SourceFiles
    kOmegaSSTLM.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTLM_H
#define kOmegaSSTLM_H

#include "RAS/kOmegaSST/kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                         Class kOmegaSSTLM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTLM
:
    public kOmegaSST<BasicTurbulenceModel>
{
protected:

    // Protected data

        // Model constants

            dimensionedScalar ca1_;
            dimensionedScalar ca2_;

            dimensionedScalar ce1_;
            dimensionedScalar ce2_;

            dimensionedScalar cThetat_;
            dimensionedScalar sigmaThetat_;
            dimensionedScalar cCross_;

            //- Convergence criterion for the lambda/thetat loop
            scalar lambdaErr_;

            //- Maximum number of iterations to converge the lambda/thetat loop
            label maxLambdaIter_;

            //- Stabilization for division by the magnitude of the velocity
            const dimensionedScalar deltaU_;

            //- Stationary crossflow correction
            Switch crossflow_;

            //- Surface roughness height (rms)
            scalar h_;

        // Fields

            //- Transition onset momentum-thickness Reynolds number
            volScalarField ReThetat_;

            //- Intermittency
            volScalarField gammaInt_;

            //- Effective intermittency
            volScalarField::Internal gammaIntEff_;


    // Protected Member Functions

        //- Modified form of the k-omega SST F1 function
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G,
            const bool deferred = false
        ) const;

        //- Modified form of the k-omega SST epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU,
            const bool deferred = false
        ) const;

        //- Freestream blending-function
        tmp<volScalarField::Internal> Fthetat
        (
            const volScalarField::Internal& Us,
            const volScalarField::Internal& Omega,
            const volScalarField::Internal& nu
        ) const;

        //- Freestream blending-function for the crossflow sink term
        tmp<volScalarField::Internal> Fthetat2
        (
            const volScalarField::Internal& Us,
            const volScalarField::Internal& Omega,
            const volScalarField::Internal& nu
        ) const;

        //- Empirical correlation for stationary crossflow
        tmp<volScalarField::Internal> ReScf
        (
            const volVectorField& Urel,
            const volScalarField::Internal& Us,
            const volScalarField::Internal& nu
        ) const;

        //- Empirical correlation for critical Reynolds number where the
        //  intermittency first starts to increase in the boundary layer
        tmp<volScalarField::Internal> ReThetac() const;

        //- Empirical correlation that controls the length of the
        //  transition region
        tmp<volScalarField::Internal> Flength
        (
            const volScalarField::Internal& nu
        ) const;

        //- Transition onset location control function
        tmp<volScalarField::Internal> Fonset
        (
            const volScalarField::Internal& Rev,
            const volScalarField::Internal& ReThetac,
            const volScalarField::Internal& RT
        ) const;

        //- Return the transition onset momentum-thickness Reynolds number
        // (based on freestream conditions)
        tmp<volScalarField::Internal> ReThetat0
        (
            const volScalarField::Internal& Us,
            const volScalarField::Internal& dUsds,
            const volScalarField::Internal& nu
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void coupledReThetatGammaIntCorrect();
        void segregatedReThetatGammaIntCorrect();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTLM");


    // Constructors

        //- Construct from components
        kOmegaSSTLM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        kOmegaSSTLM(const kOmegaSSTLM&) = delete;


    //- Destructor
    virtual ~kOmegaSSTLM()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function transition onset momentum-thickness Reynolds number
        const volScalarField& ReThetat() const
        {
            return ReThetat_;
        }

        //- Access function to intermittency
        const volScalarField& gammaInt() const
        {
            return gammaInt_;
        }

        //- Return the effective diffusivity for transition onset
        //  momentum-thickness Reynolds number
        tmp<volScalarField> DReThetatEff() const
        {
            return volScalarField::New
            (
                "DReThetatEff",
                sigmaThetat_*(this->nut_ + this->nu())
            );
        }

        //- Return the effective diffusivity for intermittency
        tmp<volScalarField> DgammaIntEff() const
        {
            return volScalarField::New
            (
                "DgammaIntEff",
                this->nut_ + this->nu()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kOmegaSSTLM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTLM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
