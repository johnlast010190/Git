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

Class
    Foam::solidThermo

Description
    Base-class for solid thermodynamic properties.

SourceFiles
    solidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermo_H
#define solidThermo_H

#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMesh/fvMesh.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "basicThermo/basicThermo.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermo
:
    virtual public basicThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("solidThermo");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        dictionary,
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName
        ),
        (obr, dict, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<solidThermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Selection using a specified dictionary
        static autoPtr<solidThermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidThermo();


    // Member functions

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const = 0;

        //- Density [kg/m^3]
        virtual tmp<volScalarField> rho() const = 0;

        //- Density for patch [kg/m^3]
        virtual tmp<scalarField> rho(const label patchi) const = 0;

        //- Density [kg/m^3] (not applicable for solids)
        virtual tmp<volScalarField> buoyantRho() const = 0;

        //- Return non-const access to the local density field [kg/m^3]
        virtual volScalarField& rho() = 0;

        //- Thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const = 0;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa(const label patchi) const = 0;

        //- Anisotropic thermal conductivity
        //  in the local coordinate system [W/m/K]
        virtual tmp<volSymmTensorField> KappaLocal() const = 0;

        //- Anisotropic thermal conductivity for patch
        //  in the local coordinate system [W/m/K]
        virtual tmp<symmTensorField> KappaLocal(const label patchi) const = 0;

        //- Anisotropic thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volSymmTensorField> alphahLocal() const = 0;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const = 0;
};


/*---------------------------------------------------------------------------*\
                Class solidThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidThermo::implementation
:
    virtual public solidThermo
{
protected:

    // Protected data

        //- Pressure [Pa]
        //  Note: This value should never be used. Solid thermo should only be
        //  instantiated with thermo models that do not depend on pressure.
        volScalarField p_;

        //- Density field [kg/m^3]
        volScalarField rho_;


public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry& obr, const word& phaseName);

        //- Construct from mesh, dictionary and phase name
        implementation
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }

        //- Pressure [Pa]
        //  Non-const access allowed for transport equations
        virtual volScalarField& p();

        //- Pressure [Pa]
        virtual const volScalarField& p() const;

        //- Density [kg/m^3]
        virtual tmp<volScalarField> rho() const;

        //- Density for patch [kg/m^3]
        virtual tmp<scalarField> rho(const label patchi) const;

        //- Return non-const access to the local density field [kg/m^3]
        virtual volScalarField& rho();

        //- Density [kg/m^3] (not applicable for solids)
        virtual tmp<volScalarField> buoyantRho() const
        {
            NotImplemented;
        }

        //- Anisotropic thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volSymmTensorField> alphahLocal() const
        {
            NotImplemented;
        }
};


/*---------------------------------------------------------------------------*\
                   Class solidThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidThermo::composite
:
    public basicThermo::implementation,
    public solidThermo::implementation
{
public:

    // Constructors

        //- Construct from mesh and phase name
        composite
        (
            const objectRegistry& obr,
            const word& phaseName
        )
        :
            basicThermo::implementation(obr, phaseName),
            solidThermo::implementation(obr, phaseName)
        {}

        //- Construct from mesh, dictionary and phase name
        composite
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName
        )
        :
            basicThermo::implementation(obr, dict, phaseName),
            solidThermo::implementation(obr, dict, phaseName)
        {}


    // Member Functions

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
