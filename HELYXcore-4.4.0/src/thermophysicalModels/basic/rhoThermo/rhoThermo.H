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
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::rhoThermo

Group
    grpRhoThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    rhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermo_H
#define rhoThermo_H

#include "fluidThermo/fluidThermo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class rhoThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;

    //- Runtime type information
    TypeName("rhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    //- Selector
    static autoPtr<rhoThermo> New
    (
        const objectRegistry&,
        const word& phaseName = word::null
    );


    //- Destructor
    virtual ~rhoThermo();


    // Member functions

        //- Update properties
        virtual void correct() = 0;

        // Derived thermodynamic properties

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho() = 0;

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
};


/*---------------------------------------------------------------------------*\
                 Class rhoThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo::implementation
:
    virtual public rhoThermo
{
protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'thermo:rho' to avoid (potential) conflict with solver density
        volScalarField rho_;


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

        //- Update properties
        virtual void correct()
        {}


        // Derived thermodynamic properties

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            //  Limit thermo rho between rhoMin and rhoMax
            virtual void correctRho
            (
                const volScalarField& deltaRho,
                const dimensionedScalar& rhoMin,
                const dimensionedScalar& rhoMax
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


/*---------------------------------------------------------------------------*\
                    Class rhoThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public rhoThermo::implementation
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
            fluidThermo::implementation(obr, phaseName),
            rhoThermo::implementation(obr, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
