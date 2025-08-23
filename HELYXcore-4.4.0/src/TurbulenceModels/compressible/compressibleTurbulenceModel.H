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
    (c) 2010-2024 Engys Ltd.
    (c) 2013-2015 OpenFOAM Foundation

Class
    Foam::compressibleTurbulenceModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    compressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTurbulenceModel_H
#define compressibleTurbulenceModel_H

#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                Class compressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleTurbulenceModel
:
    public turbulenceModel
{
protected:

    // Protected data

        const volScalarField& rho_;


    // Protected member functions

        //- ***HGW Temporary function to be removed when the run-time selectable
        //  thermal transport layer is complete
        virtual void correctNut()
        {}


public:

    //- Runtime type information
    TypeName("compressibleTurbulenceModel");


    // Constructors

        //- Construct from components
        compressibleTurbulenceModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const word& propertiesName
        );

        //- Disallow default bitwise copy construct
        compressibleTurbulenceModel
        (
            const compressibleTurbulenceModel&
        ) = delete;


    //- Destructor
    virtual ~compressibleTurbulenceModel()
    {}


    // Member functions

        //- Return the density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            volVectorField& U,
            const word& scheme
        ) const = 0;

        //- Correct the turbulence thermal diffusivity for energy transport
        virtual void correctEnergyTransport();

        // Access functions overridden by Diffusivity
        // (ThermalDiffusiviy, EddyDiffusivity, ...)

            //- Return the thermal diffusivity for temperature [J/m/s/K]
            virtual tmp<volScalarField> kappa() const
            {
                NotImplemented;
            }

            //- Return the laminar thermal diffusivity for temperature on patch
            //  [J/m/s/K]
            virtual tmp<scalarField> kappa(const label patchi) const
            {
                NotImplemented;
            }

            //- Return the laminar thermal diffusivity for enthalpy [kg/m/s]
            virtual tmp<volScalarField> alpha() const
            {
                NotImplemented;
            }

            //- Return the laminar thermal diffusivity for enthalpy on patch
            //  [kg/m/s]
            virtual tmp<scalarField> alpha(const label patchi) const
            {
                NotImplemented;
            }

            //- Return the turbulent thermal diffusivity for enthalpy [kg/m/s]
            virtual tmp<volScalarField> alphat() const
            {
                NotImplemented;
            }

            //- Return the turbulent thermal diffusivity for enthalpy for a patch
            //  [kg/m/s]
            virtual tmp<scalarField> alphat(const label patchi) const
            {
                NotImplemented;
            }

            //- Return the effective turbulent thermal diffusivity for temperature
            //  [J/m/s/K]
            virtual tmp<volScalarField> kappaEff() const
            {
                NotImplemented;
            }

            //- Return the effective turbulent thermal diffusivity for temperature
            //  [J/m/s/K]
            virtual tmp<scalarField> kappaEff(const label patchi) const
            {
                NotImplemented;
            }

            //- Return the effective turbulent thermal diffusivity for enthalpy
            //  [kg/m/s]
            virtual tmp<volScalarField> alphaEff() const
            {
                NotImplemented;
            }

            //- Return the effective turbulent thermal diffusivity for enthalpy
            //  for a patch [kg/m/s]
            virtual tmp<scalarField> alphaEff(const label patchi) const
            {
                NotImplemented;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleTurbulenceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
