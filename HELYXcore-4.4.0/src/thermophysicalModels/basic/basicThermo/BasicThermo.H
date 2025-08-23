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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021-2025 Engys Ltd.

Class
    Foam::BasicThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    BasicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef BasicThermo_H
#define BasicThermo_H

#include "materialModels/materialTables/materialTables.H"
#include "materialsSpecie/specieModels/specieModels.H"
#include "equationOfState/equationOfState/equationOfState.H"
#include "thermo/thermo/thermodynamics/thermodynamics.H"
#include "thermo/energyConversion/energyConversion/energyConversion.H"
#include "thermo/energy/energy/energy.H"
#include "thermo/energy/energyInverseTemperature/energyInverseTemperature.H"
#include "thermo/energyConversion/standardThermo/standardThermo.H"
#include "materialsSpecie/materialsSpecie.H"
#include "include/dummyThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class BasicThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermoType, class MixtureType>
class BasicThermo
:
    public BasicThermoType,
    public MixtureType
{
protected:

    // Protected data

        //- Reference to material models
        materialTables& materials_;

        //- Energy
        volScalarField he_;

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        // Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;

        //- Is the fluid buoyant?
        bool isBuoyant_;

        //- Is boyancy distinct?
        bool isDistinctBuoyancy_;

        //- Do we calcualte T from he or vice versa?
        bool calculateTFromhe_;


    // Protected Member Functions

        //- Correct the enthalpy/internal energy field boundary
        void heBoundaryCorrection(fvPatchScalarField& pf);

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);

        //- Construct material tables
        materialTables& matLookupOrConstruct
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


private:

    // Private Member Functions

        //- Load liquid specific models
        void loadLiquidModels();

        //- Initialize
        void init();

        //- Initialize all the he times
        void initHe
        (
            const volScalarField& p,
            volScalarField& he
        );

        //- Update T, p
        tmp<scalarField> tmpBoundaryField
        (
            const scalarField& T,
            const label patchi,
            const word& modelName
        ) const;

        //- Create name of the energy field based on the model
        word energyName(const dictionary& dict) const;

        //- Initialize buoyancy flags
        void initBuoyancyFlags(const word& buoyantName);

        //- Report solution flags
        void reportFlags() const;


public:

    // Constructors

        //- Construct from object registry, dictionary and phase name
        BasicThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Construct from object registry and phase name
        BasicThermo
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Return material tables
        virtual materialTables& materials() const
        {
            return materials_;
        }


    //- Destructor
    virtual ~BasicThermo();


    // Member functions

        //- Addding phasic variable to list in material properties
        //  (advantage is that it is more flexible and can be used for
        //  different solver combinations)
        virtual void addPhasicVariable(const word& var)
        {
            materials_.addPhasicVariable(var);
        }

        //- Checking if variable is phasic
        //  (advantage is that it is more flexible and can be used for
        //  different solver combinations)
        virtual bool isPhasicVariable(const word& var) const
        {
            return materials_.isPhasicVariable(var);
        }

        //- Adding phase name to variable if it is found in the list inside
        //  material properties
        virtual word phasePropertyName
        (
            const word& name,
            const word& phaseName
        ) const
        {
            return materials_.phasePropertyName(name, phaseName);
        }

        virtual word phasePropertyName(const word& name) const
        {
            return basicThermo::phasePropertyName(name);
        }

        //- Return the compostion of the mixture
        virtual typename MixtureType::basicMixtureType& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        const virtual typename MixtureType::basicMixtureType&
        composition() const
        {
            return *this;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return dummyThermo::typeName();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;

        //- Return true if Cp or Cv (depending selected energy type) is constant
        virtual bool isCpvConst() const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const;

        //- Check if all material models are const
        virtual bool isConst() const;

        //- Is the model buoyant?
        virtual bool buoyant() const
        {
            return isBuoyant_;
        }

        //- It the buoyant rho same as main rho?
        virtual bool distinctBuoyancy() const
        {
            return isDistinctBuoyancy_;
        }

        //- Do we solve for T from he or vice-versa?
        virtual bool calculatesTFromhe() const
        {
            return calculateTFromhe_;
        }

        //- Set whether to solve for T from he or vice versa
        virtual void setToCalculateTFromhe(bool calculateTFromhe)
        {
            calculateTFromhe_ = calculateTFromhe;
        }


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const
            {
                return Cp_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const
            {
                return Cv_;
            }

            //- Enthalpy
            bool isEnthalpy() const;


        // Derived thermodynamic properties

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

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy for patch [J/kg/K]
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
            ) const
            {
                // Needed only for legacy thermo
                NotImplemented;
            }

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& he,
                const scalarField& T0,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const;


        // Fields derived from transport state variables

            //- Effective thermal turbulent conductivity of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

            //- Effective thermal turbulent conductivity of mixture
            //  for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective turbulent thermal diffusivity of energy
            //  of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective turbulent thermal diffusivity of energy
            //  of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Density [kg/m^3] (density for use in buoyancy term)
            virtual tmp<volScalarField> buoyantRho() const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "BasicThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
