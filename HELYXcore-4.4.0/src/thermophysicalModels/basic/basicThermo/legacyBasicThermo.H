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

Class
    Foam::legacyBasicThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    legacyBasicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef legacyBasicThermo_H
#define legacyBasicThermo_H

#include "mixtures/basicMixture/basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class legacyBasicThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class MixtureType>
class legacyBasicThermo
:
    public BasicThermo,
    public MixtureType
{
protected:

    // Protected data

        //- Energy field
        volScalarField he_;

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        // Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;


    // Protected Member Functions

        //- Return a volScalarField of the given property
        template
        <
            class CellMixture,
            class PatchFaceMixture,
            class Method,
            class ... Args
        >
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            CellMixture cellMixture,
            PatchFaceMixture patchFaceMixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class PatchFaceMixture, class Method, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            PatchFaceMixture patchFaceMixture,
            Method psiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Initialize legacyBasicThermo
        void init
        (
            const volScalarField& p,
            const volScalarField& T,
            volScalarField& he
        );

        //- Correct the enthalpy/internal energy field boundary
        void heBoundaryCorrection(fvPatchScalarField& pf);

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);


public:

    // Constructors

        //- Construct from obr
        legacyBasicThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from obr and dictionary
        legacyBasicThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        legacyBasicThermo
        (
            const legacyBasicThermo<BasicThermo, MixtureType>&
        ) = delete;


    //- Destructor
    virtual ~legacyBasicThermo();


    // Member Functions

        //- Return the composition of the mixture
        virtual typename MixtureType::basicMixtureType& composition()
        {
            return *this;
        }

        //- Return the composition of the mixture
        virtual const typename MixtureType::basicMixtureType&
        composition() const
        {
            return *this;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return MixtureType::thermoType::typeName();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return MixtureType::thermoType::incompressible;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return MixtureType::thermoType::isochoric;
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
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
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

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "basicThermo/legacyBasicThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
