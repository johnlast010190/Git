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
    (c) 2025 Engys Ltd.

Class
    Foam::RedlichKwongPengRobinsonGas

Description
    RedlichKwongPengRobinsonGas cubic equation of state for gases
    Coefficient mixing is very inaccurate and not supported, so this equation of
    state can only be applied to mixtures by using a thermophysical properties
    model which mixes specie property value.

    Source: Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M.,
    Thermodynamic modeling for numerical simulations based on
    the generalized cubic equation of state.

Usage
    \table
        Property     | Description
        Tc           | Critical temperature [K]
        Vc           | Critical volume [m^3/kmol]
        Pc           | Critical pressure [Pa]
        omega        | Acentric factor [-]
    \endtable

    Example specification of the RedlichKwongPengRobinsonGas
    equation of state for methane:

    \verbatim
    equationOfState RedlichKwongPengRobinsonGas;
    equationOfStateCoeffs
    {
        Tc        190.55;
        Vc        0.0285;
        Pc        4.595e6;
        omega     0.0115;
    }
    \endverbatim

SourceFiles
    RedlichKwongPengRobinsonGasI.H
    RedlichKwongPengRobinsonGas.C

\*---------------------------------------------------------------------------*/

#ifndef RedlichKwongPengRobinsonGas_H
#define RedlichKwongPengRobinsonGas_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RedlichKwongPengRobinsonGas Declaration
\*---------------------------------------------------------------------------*/

class RedlichKwongPengRobinsonGas
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

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;

        //- kappa /c_alpha for a*alpha term in EoS [-]
        scalar kappa_;

        //- a coefficient and A for EoS [-]
        scalar a_;

        //- b coefficient for Eos [-]
        scalar b_;

        // - delta1 (d1), delta2 (d2) and short cuts
        // for generalized EdeoS coeff.[-]
        scalar d1_, d2_, d1PlusD2_, d1TimesD2_, d1MinusD2_;

        //- Quick model access
        enum {Rind, Wind, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary and object registry
        RedlichKwongPengRobinsonGas
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<RedlichKwongPengRobinsonGas> clone() const;

        //- Selector from dictionary
        inline static autoPtr<RedlichKwongPengRobinsonGas> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("RedlichKwongPengRobinsonGas");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override;

        //- Return alpha
        inline scalar alpha(scalar T) const;

        //- Return d(a*alpha)/dT
        inline scalar ap(scalar T) const;

        //- Return d2(a*alpha)/dT2
        inline scalar app(scalar T) const;

        //- Return A
        inline scalar A(scalar p, scalar T) const;

        //- Return B
        inline scalar B(scalar p, scalar T) const;

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsThree(rho, p_, T_, sMod_(Rind));

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsThree(hContribution, p_, T_, sMod_(Rind));

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsThree(CpContribution, p_, T_, sMod_(Wind));

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsThree(eContribution, p_, T_, sMod_(Rind));

        //- Functions for Cv contribution [J/kg/K]
        scalarFieldFunctionsThree(CvContribution, p_, T_, sMod_(Wind));

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsThree(spContribution, p_, T_, sMod_(Rind));

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsZero(svContribution);

        //- Functions for compressibility [s^2/m^2]
        scalarFieldFunctionsThree(psi, p_, T_, sMod_(Rind));

        //- Functions for compression factor []
        scalarFieldFunctionsTwo(Z, p_, T_);

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsThree(CpMCv, p_, T_, sMod_(Rind));

        //- Functions for volumetric coefficient of thermal expansion [1/T]
        scalarFieldFunctionsZero(alphav);

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

#include "RedlichKwongPengRobinsonGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
