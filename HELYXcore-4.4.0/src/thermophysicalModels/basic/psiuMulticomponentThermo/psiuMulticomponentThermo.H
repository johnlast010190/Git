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
    (c) 2025 Engys Ltd.

Class
    Foam::psiuMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    psiuMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiuMulticomponentThermo_H
#define psiuMulticomponentThermo_H

#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "mixtures/basicCombustionMixture/basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiuMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo
:
    virtual public rhoThermo
{
protected:

    // Protected Member Functions

        wordList heuBoundaryTypes();
        void heuBoundaryCorrection(volScalarField& heu);


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("psiuMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            psiuMulticomponentThermo,
            objectRegistry,
            (const objectRegistry& obr, const word& phaseName),
            (obr, phaseName)
        );


    // Selectors

        //- Standard selection based on objectRegistry
        static autoPtr<psiuMulticomponentThermo> New
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~psiuMulticomponentThermo();


    // Member functions

        //- Return the composition of the combustion mixture
        virtual basicCombustionMixture& composition() = 0;

        //- Return the composition of the combustion mixture
        virtual const basicCombustionMixture& composition() const = 0;


        // Access to thermodynamic state variables.

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu() = 0;

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const = 0;


        // Derived thermodynamic properties

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const = 0;

            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const = 0;

            //- Unburnt gas density [kg/m^3]
            virtual tmp<volScalarField> rhou() const
            {
                return p()*psiu();
            }

            //- Burnt gas density [kg/m^3]
            virtual tmp<volScalarField> rhob() const
            {
                return p()*psib();
            }

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const = 0;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const = 0;

            //- Dynamic viscosity of unburnt gas [kg/ms]
            virtual tmp<volScalarField> muu() const = 0;

            //- Dynamic viscosity of burnt gas [kg/ms]
            virtual tmp<volScalarField> mub() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class psiuMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo::implementation
:
    virtual public psiuMulticomponentThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const objectRegistry&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                Class psiuMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public rhoThermo::implementation,
    public psiuMulticomponentThermo::implementation
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
            rhoThermo::implementation(obr, phaseName),
            psiuMulticomponentThermo::implementation(obr, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
