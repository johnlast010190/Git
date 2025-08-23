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
    (c) 2025 Engys Ltd.

Class
    Foam::PsiuMulticomponentThermo

Description
    Thermo implementation based on compressibility with additional unburnt
    thermodynamic state

SourceFiles
    PsiuMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef PsiuMulticomponentThermo_H
#define PsiuMulticomponentThermo_H

#include "basicThermo/BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PsiuMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class PsiuMulticomponentThermo
:
    public BasicThermo<BasicPsiThermo, MixtureType>
{
    // Private Data

        //- Unburnt temperature
        volScalarField Tu_;

        //- Unburnt energy
        volScalarField heu_;

        //- Types of mixtures
        enum {reactants, products};

        //- Mixture type
        const word mixtureModel_;


    // Private Member Functions

        //- Material model
        tmp<scalarField> material
        (
            const word& name,
            label mixtureType,
            label patchi = -1
        ) const;

        //- Temperature from unburned he
        tmp<scalarField> boundaryThe
        (
            const scalarField& he,
            label mixtureType,
            label patchi
        ) const;

        //- Temperature from unburned he
        tmp<scalarField> internalThe
        (
            const scalarField& he,
            label mixtureType
        ) const;

        //- Calculate the thermo variables
        void calculate();

        //- Materials object registry
        const objectRegistry& matDb() const
        {
            return this->T_.db().subRegistry("materialModels");
        }

        //- Mixture fraction
        const volScalarField& ft()
        {
            return this->T_.db().template lookupObject<volScalarField>("ft");
        }

        //- Regress fraction
        const volScalarField& b()
        {
            return this->T_.db().template lookupObject<volScalarField>("b");
        }


public:

        //- Runtime type information
        TypeName("premixed");

    // Constructors

        //- Construct from mesh and phase name
        PsiuMulticomponentThermo
        (
            const objectRegistry& obr,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        PsiuMulticomponentThermo
        (
            const PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>&
        ) = delete;


    //- Destructor
    virtual ~PsiuMulticomponentThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // Thermodynamic state

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const
            {
                return Tu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const
            {
                return heu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu()
            {
                return heu_;
            }


        // Derived thermodynamic properties

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const;


        // Derived transport properties

            //- Dynamic viscosity of unburnt gas [kg/m/s]
            virtual tmp<volScalarField> muu() const;

            //- Dynamic viscosity of burnt gas [kg/m/s]
            virtual tmp<volScalarField> mub() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const PsiuMulticomponentThermo<BasicPsiThermo, MixtureType>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "PsiuMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
