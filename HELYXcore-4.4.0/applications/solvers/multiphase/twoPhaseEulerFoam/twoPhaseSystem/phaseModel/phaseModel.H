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
    (c) 2010-2016 Engys Ltd.
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "db/dictionary/dictionary.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "transportModel/transportModel.H"
#include "rhoThermo/rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class twoPhaseSystem;
class diameterModel;

template<class Phase>
class PhaseCompressibleTurbulenceModel;


/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField,
    public transportModel
{
    // Private data

        //- Reference to the twoPhaseSystem to which this phase belongs
        const twoPhaseSystem& fluid_;

        //- Name of phase
        word name_;

        dictionary phaseDict_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Thermophysical properties
        autoPtr<rhoThermo> thermo_;

        //- Velocity
        volVectorField U_;

        //- Volumetric flux of the phase
        surfaceScalarField alphaPhi_;

        //- Mass flux of the phase
        surfaceScalarField alphaRhoPhi_;

        //- Volumetric flux of the phase
        autoPtr<surfaceScalarField> phiPtr_;

        //- Diameter model
        autoPtr<diameterModel> dPtr_;

        //- Turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel>> turbulence_;


public:

    // Constructors

        phaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the twoPhaseSystem to which this phase belongs
        const twoPhaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Return the other phase in this two-phase system
        const phaseModel& otherPhase() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Return the turbulence model
        const PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence() const;

        //- Return non-const access to the turbulence model
        //  for correction
        PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence();

        //- Return the thermophysical model
        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        //- Return non-const access to the thermophysical model
        //  for correction
        rhoThermo& thermo()
        {
            return thermo_();
        }

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return thermo_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return thermo_->nu(patchi);
        }

        //- Return the laminar dynamic viscosity
        const volScalarField& mu() const
        {
            return thermo_->mu();
        }

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        const volScalarField& kappa() const
        {
            return thermo_->kappa();
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->kappaEff(alphat);
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->kappaEff(alphat, patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->alphaEff(alphat);
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->alphaEff(alphat, patchi);
        }

        //- Return the specific heat capacity
        const volScalarField& Cp() const
        {
            return thermo_->Cp();
        }

        //- Heat capacity at constant pressure/volume [J/kg/K]
        const volScalarField& Cpv() const
        {
            return thermo_->Cpv();
        }

        //- Return the density
        tmp<volScalarField> rho() const
        {
            return thermo_->rho();
        }

        //- dummy return conductivity
        tmp<volScalarField> lambda() const
        {
            NotImplemented;
        }

        //- dummy return Prt
        tmp<volScalarField> Prt() const
        {
            NotImplemented;
        }

        //- dummy return Pr
        tmp<volScalarField> Pr() const
        {
            NotImplemented;
        }

        //- dummy return alphaLam
        tmp<volScalarField> alphaLam() const
        {
            NotImplemented;
        }

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return non-const access to the velocity
        //  Used in the momentum equation
        volVectorField& U()
        {
            return U_;
        }

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Return non-const access to the volumetric flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return the volumetric flux of the phase
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Return non-const access to the volumetric flux of the phase
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Return the mass flux of the phase
        const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return non-const access to the mass flux of the phase
        surfaceScalarField& alphaRhoPhi()
        {
            return alphaRhoPhi_;
        }

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        //- Correct the phase properties
        //  other than the thermodynamics and turbulence
        //  which have special treatment
        void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

        //- Dummy Read for transportModel
        virtual bool read()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
