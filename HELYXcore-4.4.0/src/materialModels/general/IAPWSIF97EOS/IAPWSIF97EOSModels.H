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
    (c) 2022-2025 Engys Ltd.

Class
    Foam::IAPWSIF97EOSModels

Description
    This calculation page is based on the Thermodynamic Properties of Water
    and Steam (Reference) [1] (named as IAPWS-IF97) provided by IAPWS.
    Detailed information about used equations, constants, range of validity
    etc is presented in PDF version of IAPWS.

    Release which can be downloaded from IAPWS web site www.iapws.org

Note:
    1. These models are not fully consistent with the material library.
    Additional work is required to fully integrate.
    Now the gas constant R is hard-coded + in case of computing Cp
    the rho model is requred for region3 that is used directly from
    standard.
    2. No contribution models used.

Reference
    1. Revised Release on the IAPWS Industrial Formulation 1997 for
    the Thermodynamic Properties of Water and Steam.
    Lucerne, Switzerland, August 2007, available at http://www.iapws.org
    2. Kretzschmar H, Wagner W International Steam Tables.

SourceFiles
    IAPWSIF97EOSModelsI.H
    IAPWSIF97EOSModels.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97EOSModels_H
#define IAPWSIF97EOSModels_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"
#include "general/IAPWSIF97EOS/IAPWSIF97EOS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IAPWSIF97EOSModels Declaration
\*---------------------------------------------------------------------------*/

class IAPWSIF97EOSModels
:
    public materialModel,
    public rhoModel,
    public psiModel,
    public ZModel,
    public CpMCvModel,
    public CpModel,
    public CvModel,
    public haModel,
    public hsModel,
    public hfModel,
    public sModel,
    public muModel,
    public kappaModel,
    public CpContributionModel,
    public hContributionModel,
    public spContributionModel,
    public CvContributionModel,
    public eContributionModel,
    public svContributionModel
{
    // Private member data

        //- Model
        IAPWSIF97EOS waterSteamModel_;

        //- Enthalpy of formation [J/kg]
        scalar hf_;

        //- Entropy [J/kg/K]
        scalar Sf_;

        //- Should the model be allowed to cross to saturation curve?
        bool supportPhaseChange_;

        //- Region number
        label regionNumber_;

        //- Checking if we are not crossing to other regions
        inline void checkPhaseChange
        (
            const scalar p,
            const scalar T
        ) const;


public:

    // Constructors

        //- Construct from dictionary and object registry
        IAPWSIF97EOSModels
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<IAPWSIF97EOSModels> clone() const;

        //- Selector from dictionary
        inline static autoPtr<IAPWSIF97EOSModels> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("IAPWSIF97EOS");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return false;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return false;
        }

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsTwo(rho, p_, T_);

        //- Functions for compressibility [s^2/m^2]
        scalarFieldFunctionsTwo(psi, p_, T_);

        //- Functions for compression factor []
        scalarFieldFunctionsNone(Z);

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsTwo(CpMCv, p_, T_);

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsTwo(ha, p_, T_);

        //- Sensible enthalpy [J/kg]
        scalarFieldFunctionsTwo(hs, p_, T_);

        //- Enthalpy of formation [J/kg]
        scalarFieldFunctionsZero(hf);

        //- Entropy [J/kg/K]
        scalarFieldFunctionsTwo(s, p_, T_);

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsTwo(Cp, p_, T_);

        //- Functions for Cv [J/kg/K]
        scalarFieldFunctionsTwo(Cv, p_, T_);

        //- Dynamic viscosity [kg/ms]
        scalarFieldFunctionsTwo(mu, p_, T_);

        //- Thermal conductivity [W/mK]
        scalarFieldFunctionsTwo(kappa, p_, T_);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsNone(CpContribution);

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsNone(hContribution);

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsNone(spContribution);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsNone(CvContribution);

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsNone(eContribution);

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsNone(svContribution);

        //- Return model dependence
        materialDepList();

        //- Update dependent model table
        virtual void updateTable(const word&) override
        {}

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

#include "IAPWSIF97EOSModelsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
