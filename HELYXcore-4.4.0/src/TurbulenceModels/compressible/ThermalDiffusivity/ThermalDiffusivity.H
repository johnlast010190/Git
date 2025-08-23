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
    (c) 2015-2022 OpenFOAM Foundation

Class
    Foam::ThermalDiffusivity

Description
    Templated wrapper class to provide compressible turbulence models
    thermal diffusivity based thermal transport.

SourceFiles
    ThermalDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalDiffusivity_H
#define ThermalDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class ThermalDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class ThermalDiffusivity
:
    public BasicTurbulenceModel
{

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef volScalarField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct
        ThermalDiffusivity
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<ThermalDiffusivity> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


        //- Return a reference to the selected turbulence model
        static autoPtr<ThermalDiffusivity> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~ThermalDiffusivity()
    {}


    // Member Functions

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappa() const
        {
            return volScalarField::New("kappa", this->transport_.kappa());
        }

        //- Thermal diffusivity for temperature of mixture
        //  for patch [J/m/s/K]
        virtual tmp<scalarField> kappa(const label patchi) const
        {
            return this->transport_.kappa().boundaryField()[patchi];
        }

        //- Turbulent thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alphat() const;

        //- Turbulent thermal diffusivity for enthalpy for a patch [kg/m/s]
        virtual tmp<scalarField> alphat(const label patchi) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return kappa();
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return kappa(patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return volScalarField::New
            (
                "alphaEff",
                this->transport_.kappa()/this->transport_.Cpv()
            );
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return
                this->transport_.kappa().boundaryField()[patchi]
               /this->transport_.Cpv().boundaryField()[patchi];
        }

        //- Return the turbulent mass diffusivity [kg/m/s]
        virtual tmp<volScalarField> Dmt() const;

        //- Return the turbulent mass diffusivity for a patch [kg/m/s]
        virtual tmp<scalarField> Dmt(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermalDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
