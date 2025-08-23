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
    (c) 2012-2022 OpenFOAM Foundation

Class
    Foam::fluidThermo

Description
    Fundamental fluid thermodynamic properties

SourceFiles
    fluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    virtual public basicThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const objectRegistry&,
        const word& phaseName = word::null
    );


    //- Destructor
    virtual ~fluidThermo();


    // Member functions

        // Access to thermodynamic state variables

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            //  Limit thermo rho between rhoMin and rhoMax
            virtual void correctRho
            (
                const volScalarField& deltaRho,
                const dimensionedScalar& rhoMin,
                const dimensionedScalar& rhoMax
            ) = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const = 0;

            //- Gauge pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p() = 0;

            //- Gauge pressure [Pa]
            virtual const volScalarField& p() const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;


        // Fields derived from transport state variables

            //- Effective thermal turbulent conductivity of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal turbulent conductivity of mixture
            //  for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective turbulent thermal diffusivity of energy
            //  of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective turbulent thermal diffusivity of energy
            //  of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Density [kg/m^3]
            // (rho for buoyancy therm, available only in material library)
            virtual tmp<volScalarField> buoyantRho() const = 0;

            virtual bool isConst() const = 0;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;
};


/*---------------------------------------------------------------------------*\
                Class fluidThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo::implementation
:
    virtual public fluidThermo
{
protected:

    // Protected Data

        //- Pressure [Pa]
        volScalarField& p_;

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;


public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);

        //- Construct from mesh, dictionary and phase name
        implementation
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Pressure [Pa]
        //  Non-const access allowed for transport equations
        virtual volScalarField& p();

        //- Pressure [Pa]
        virtual const volScalarField& p() const;

        //- Compressibility [s^2/m^2]
        virtual const volScalarField& psi() const;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual const volScalarField& mu() const;

        //- Density [kg/m^3]
        // (rho for buoyancy therm, available only in material library)
        virtual tmp<volScalarField> buoyantRho() const
        {
            NotImplemented;
        }

        virtual bool isConst() const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
