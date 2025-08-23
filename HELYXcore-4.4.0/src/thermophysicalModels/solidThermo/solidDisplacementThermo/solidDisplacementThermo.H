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
    (c) 2019-2023 OpenFOAM Foundation

Class
    Foam::solidDisplacementThermo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidDisplacementThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacementThermo_H
#define solidDisplacementThermo_H

#include "solidThermo/solidThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class solidDisplacementThermo Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacementThermo
:
    virtual public solidThermo::composite
{
    // Private data

        //- Reference to material models
        materialTables& materials_;

        //- Switch to enable plane stress
        Switch planeStress_;

        //- Switch to enable thermal stress
        Switch thermalStress_;

        //- Heat capacity at constant pressure [J/kg/K]
        volScalarField Cp_;

        //- Thermal diffusivity for temperature [W/m/K]
        volScalarField kappa_;

        //- Youngs modulus [Pa]
        volScalarField E_;

        //- Poisson's ratio []
        volScalarField nu_;

        //- Volumetric thermal expansion coefficient [1/T]
        volScalarField alphav_;


    // Private member functions

        void loadModels(volScalarField& prop) const;

        materialTables& matLookupOrConstruct
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


public:

    //- Runtime type information
    TypeName("solidDisplacementThermo");


    // Constructors

        //- Construct from mesh and phase name
        solidDisplacementThermo
        (
            const objectRegistry& obr,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidDisplacementThermo();


    // Member Functions

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return type();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return true;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return true;
        }

        //- Returns true to enable plane stress
        bool planeStress() const
        {
            return planeStress_;
        }

        //- Returns true to enable thermal stress
        bool thermalStress() const
        {
            return thermalStress_;
        }


        // Access to thermophysical state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Youngs modulus [Pa]
            virtual const volScalarField& E() const;

            //- Youngs modulus for a patch [Pa]
            virtual const scalarField& E(const label patchi) const;

            //- Poisson's ratio []
            virtual const volScalarField& nu() const;

            //- Poisson's ratio for a patch[]
            virtual const scalarField& nu(const label patchi) const;

            //- Volumetric thermal expansion coefficient [1/T]
            virtual const volScalarField& alphav() const;

            //- Volumetric thermal expansion coefficient for a patch [1/T]
            virtual const scalarField& alphav(const label patchi) const;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                NotImplemented;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                NotImplemented;
            }


        // Derived thermodynamic properties

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const
            {
                NotImplemented;
            }

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const
            {
                NotImplemented;
            }

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const
            {
                NotImplemented;
            }

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const
            {
                NotImplemented;
            }

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const
            {
                NotImplemented;
            }

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const
            {
                NotImplemented;
            }

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const
            {
                NotImplemented;
            }

            //- Absolute enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hf() const
            {
                NotImplemented;
            }

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> The
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const
            {
                NotImplemented;
            }

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual const volScalarField& kappa() const;

            //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff(const volScalarField&) const
            {
                NotImplemented;
            }

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const
            {
                NotImplemented;
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                NotImplemented;
            }

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const
            {
                return true;
            }

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const
            {
                NotImplemented;
            }

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa(const label patchi) const
            {
                NotImplemented;
            }

            //- Anisotropic thermal conductivity
            //  in the local coordinate system [W/m/K]
            virtual tmp<volSymmTensorField> KappaLocal() const
            {
                NotImplemented;
            }

            //- Anisotropic thermal conductivity for patch
            //  in the local coordinate system [W/m/K]
            virtual tmp<symmTensorField> KappaLocal(const label patchi) const
            {
                NotImplemented;
            }

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p()
            {
                NotImplemented;
            }

            //- Pressure [Pa]
            virtual const volScalarField& p() const
            {
                NotImplemented;
            }


        //- Update properties
        virtual void correct();


    // I-O

        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //