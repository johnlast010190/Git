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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::liquidProperties

Description
    The thermophysical properties of a liquid.
    Generic thermophysical properties class for a liquid in which the
    functions and coefficients for each property are run-time selected.

SourceFiles
    liquidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef liquidProperties_H
#define liquidProperties_H

#include "thermophysicalProperties/thermophysicalProperties.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function2/APIdiffCoef/APIdiffCoef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class liquidProperties;
Ostream& operator<<(Ostream& os, const liquidProperties& l);

/*---------------------------------------------------------------------------*\
                    Class liquidProperties Declaration
\*---------------------------------------------------------------------------*/

class liquidProperties
:
    public thermophysicalProperties
{
    // Private data

        //- Name of the liquid
        const word name_;

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical pressure [Pa]
        scalar Pc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compressibility factor []
        scalar Zc_;

        //- Triple point temperature [K]
        scalar Tt_;

        //- Triple point pressure [Pa]
        scalar Pt_;

        //- Normal boiling temperature [K]
        scalar Tb_;

        //- Dipole moment []
        scalar dipm_;

        //- Pitzer's accentric factor []
        scalar omega_;

        //- Solubility parameter [(J/m^3)^0.5]
        scalar delta_;

        //- Rerference pressure [Pa]
        scalar pRef_;

        autoPtr<Function1<scalar>> rho_;
        autoPtr<Function1<scalar>> pv_;
        autoPtr<Function1<scalar>> hl_;
        autoPtr<Function1<scalar>> Cp_;
        autoPtr<Function1<scalar>> h_;
        autoPtr<Function1<scalar>> Cpg_;
        autoPtr<Function1<scalar>> B_;
        autoPtr<Function1<scalar>> mu_;
        autoPtr<Function1<scalar>> mug_;
        autoPtr<Function1<scalar>> kappa_;
        autoPtr<Function1<scalar>> kappag_;
        autoPtr<Function1<scalar>> sigma_;
        autoPtr<Function2s::APIdiffCoef> D_;
        autoPtr<Function1<scalar>> Dfunc1_;

        //- Liquid heat of formation [J/kg]
        scalar hf_;


public:

    TypeName("liquid");


    // Constructors

        //- Construct from dictionary
        liquidProperties(const dictionary& dict);

        //- Construct as a copy
        liquidProperties(const liquidProperties& l);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new liquidProperties(*this));
        }


    // Selectors

        //- Return a pointer to a new liquidProperties created from name
        static autoPtr<liquidProperties> New(const word& name);

        //- Return a pointer to a new liquidProperties created from dictionary
        static autoPtr<liquidProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~liquidProperties()
    {}


    // Static data

        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible = true;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric = false;


    // Member Functions

        //- Extracting fluid name from the dictionary name
        const word extractFluidName(const dictionary& dict) const;

        //- Return the name of the liquid
        virtual const word& name() const;


        // Physical constants which define the specie

            //- Mass fraction of this specie in mixture
            //  Note Mixing of liquidProperties is not currently supported
            //  so Y = 1
            inline scalar Y() const;

            //- Critical temperature [K]
            inline scalar Tc() const;

            //- Critical pressure [Pa]
            inline scalar Pc() const;

            //- Critical volume [m^3/kmol]
            inline scalar Vc() const;

            //- Critical compressibilty factor
            inline scalar Zc() const;

            //- Triple point temperature [K]
            inline scalar Tt() const;

            //- Triple point pressure [Pa]
            inline scalar Pt() const;

            //- Normal boiling temperature [K]
            inline scalar Tb() const;

            //- Dipole moment []
            inline scalar dipm() const;

            //- Pitzer's ascentric factor []
            inline scalar omega() const;

            //- Solubility parameter [(J/m^3)^(1/2)]
            inline scalar delta() const;

            //- Limit the temperature to be in the range Tlow_ to Thigh_
            inline scalar limit(const scalar T) const;


        // Fundamental equation of state properties

            //- Liquid compressibility [s^2/m^2]
            //  Note: currently it is assumed the liquid is incompressible
            inline scalar psi(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            //  Note: currently it is assumed the liquid is incompressible
            //  so CpMCv = 0
            inline scalar CpMCv(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Enthalpy contribution [J/kg]
            inline scalar hContribution(scalar p, scalar T) const;

            //- Cp contribution [J/kg/K]
            inline scalar CpContribution(scalar p, scalar T) const;

            //- Liquid sensible enthalpy [J/kg]
            inline scalar hs(scalar p, scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;

            //- Liquid absolute enthalpy [J/kg]
            inline scalar ha(scalar p, scalar T) const;

            //- Sensible internal energy [J/kg]
            //  currently pressure effects are neglected
            inline scalar es(scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            //  currently pressure effects are neglected
            inline scalar ea(scalar p, scalar T) const;

            //- Liquid entropy [J/kg/K]
            scalar s(const scalar p, const scalar T) const;


        // Physical properties

            //- Liquid density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar p, scalar T) const;

            //- Heat of vapourisation [J/kg]
            virtual scalar hl(scalar p, scalar T) const;

            //- Liquid heat capacity [J/kg/K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Vapour heat capacity [J/kg/K]
            virtual scalar Cpg(scalar p, scalar T) const;

            //- Second Virial Coefficient [m^3/kg]
            inline scalar B(scalar p, scalar T) const;

            //- Liquid viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const;

            //- Vapour viscosity [Pa s]
            virtual scalar mug(scalar p, scalar T) const;

            //- Liquid thermal conductivity  [W/(m K)]
            virtual scalar kappa(scalar p, scalar T) const;

            //- Vapour thermal conductivity  [W/(m K)]
            virtual scalar kappag(scalar p, scalar T) const;

            //- Surface tension [N/m]
            virtual scalar sigma(scalar p, scalar T) const;

            //- Vapour diffusivity [m2/s]
            virtual scalar D(scalar p, scalar T) const;

            //- Vapour diffusivity [m2/s] with specified binary pair
            virtual scalar D(scalar p, scalar T, scalar Wb) const;

            //- Invert the vapour pressure relationship to retrieve the
            //  boiling temperature as a function of pressure
            virtual scalar pvInvert(scalar p) const;


    // I-O

        //- Set value in the dictionary
        static void setDictValue
        (
            dictionary& liquidDict,
            const dictionary& dict,
            const word& key,
            const word& subDictName = word::null
        );

        //- Merge two default with user specified properties
        static dictionary mergeDicts
        (
            const dictionary& dict,
            const word& name
        );

        //- Read and set the properties present it the given dictionary
        void readIfPresent(const dictionary& dict);


        //- Write the function coefficients
        virtual void write(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const liquidProperties& l);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidProperties/liquidProperties/liquidPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
