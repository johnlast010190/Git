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
    (c) 2013-2023 OpenFOAM Foundation

Class
    Foam::twoPhaseMixtureThermo

Description

SourceFiles
    twoPhaseMixtureThermoI.H
    twoPhaseMixtureThermo.C
    twoPhaseMixtureThermoIO.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixtureThermo_H
#define twoPhaseMixtureThermo_H

#include "rhoThermo/rhoThermo.H"
#include "twoPhaseMixture/twoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class twoPhaseMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixtureThermo
:
    public rhoThermo::composite,
    public twoPhaseMixture,
    public interfaceProperties
{
    // Private data

        //- Thermo-package of phase 1
        autoPtr<rhoThermo> thermo1_;

        //- Thermo-package of phase 2
        autoPtr<rhoThermo> thermo2_;

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        // Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;


public:

    //- Runtime type information
    TypeName("twoPhaseMixtureThermo");


    // Constructors

        //- Construct from components
        twoPhaseMixtureThermo
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~twoPhaseMixtureThermo();


    // Member Functions

        const rhoThermo& thermo1() const
        {
            return thermo1_();
        }

        const rhoThermo& thermo2() const
        {
            return thermo2_();
        }

        rhoThermo& thermo1()
        {
            return thermo1_();
        }

        rhoThermo& thermo2()
        {
            return thermo2_();
        }

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Update mixture properties
        virtual void correct();

        //- Return the name of the thermo physics
        virtual word thermoName() const;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;


        // Access to thermodynamic state variables

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

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> ha() const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy for patch [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hf() const;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> The
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        // Derived thermodynamic properties

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

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual const volScalarField& kappa() const;

            //- Thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const
            {
                NotImplemented;
            }


    // IO

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
