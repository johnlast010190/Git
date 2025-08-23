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
    (c) 2015-2021 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::BoussinesqLaw

Description
    Incompressible gas equation of state using the Boussinesq approximation for
    the density as a function of temperature only:

    \verbatim
        rho = rho0*(1 - beta*(T - T0))
    \endverbatim

Usage
    \table
        Property     | Description
        rho0         | Reference density
        T0           | Reference temperature
        beta         | Coefficient of thermal expansion
    \endtable

    Example specification of the Boussinesq equation of state:

    \verbatim
    equationOfState Boussinesq;
    equationOfStateCoeffs
    {
        rho0            1;
        T0              300;
        beta            3e-03;
    }
    \endverbatim

SourceFiles
    BoussinesqLawI.H
    BoussinesqLaw.C

\*---------------------------------------------------------------------------*/

#ifndef BoussinesqLaw_H
#define BoussinesqLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BoussinesqLaw Declaration
\*---------------------------------------------------------------------------*/

class BoussinesqLaw
:
    public materialModel,
    public rhoModel,
    public hContributionModel,
    public CpContributionModel,
    public eContributionModel,
    public CvContributionModel,
    public spContributionModel,
    public svContributionModel,
    public psiModel,
    public ZModel,
    public CpMCvModel,
    public alphavModel
{
    // Private member data

        //- Reference density
        scalar rho0_;

        //- Reference temperature
        scalar T0_;

        //- Thermal expansion coefficient
        scalar beta_;

        //- Quick model access
        enum {RInd, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary and object registry
        BoussinesqLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<BoussinesqLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<BoussinesqLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("Boussinesq");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return true;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return T_->isConst();
        }

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsOne(rho, T_);

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsTwo(hContribution, p_, T_);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsTwo(CpContribution, p_, T_);

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsNone(eContribution);

        //- Functions for Cv contribution [J/kg/K]
        scalarFieldFunctionsNone(CvContribution);

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsNone(spContribution);

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsNone(svContribution);

        //- Functions for compressibility drho/dp [s^2/m^2]
        scalarFieldFunctionsNone(psi);

        //- Functions for compression factor []
        scalarFieldFunctionsThree(Z, p_, T_, sMod_(RInd));

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsNone(CpMCv);

        //- Functions for volumetric coefficient of thermal expansion [1/T]
        scalarFieldFunctionsOne(alphav, T_);

        //- Return model dependence
        materialDepList();

        //- Update dependent model table
        virtual void updateTable(const word&) override;

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoussinesqLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
