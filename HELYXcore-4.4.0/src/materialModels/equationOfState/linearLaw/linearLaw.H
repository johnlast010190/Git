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
    (c) 2013-2021 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::linearLaw

Description
    Linear equation of state with constant compressibility:

    \verbatim
        rho = rho0 + psi*p
    \endverbatim

    Coefficient mixing is very inaccurate and not supported,
    so this equation of state is not applicable to mixtures.

Usage
    \table
        Property     | Description
        rho0         | Reference density
        psi          | Constant compressibility
    \endtable

    Example specification of the linear equation of state:

    \verbatim
    equationOfState linearLaw;
    equationOfStateCoeffs
    {
        rho0        1000;
        psi         1e-5;
    }
    \endverbatim

SourceFiles
    linearLawI.H
    linearLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearLaw_H
#define linearLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class linearLaw Declaration
\*---------------------------------------------------------------------------*/

class linearLaw
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

        //- Compressibility
        scalar psi_;

        //- The reference density
        scalar rho0_;

        //- Quick model access
        enum {RInd, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary and object registry
        linearLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<linearLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<linearLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("linear");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override;

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsOne(rho, p_);

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsOne(hContribution, p_);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsNone(CpContribution);

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsNone(eContribution);

        //- Functions for Cv contribution [J/kg/K]
        scalarFieldFunctionsNone(CvContribution);

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsTwo(spContribution, p_, T_);

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsNone(svContribution);

        //- Functions for compressibility drho/dp [s^2/m^2]
        scalarFieldFunctionsZero(psi);

        //- Functions for compression factor []
        scalarFieldFunctionsThree(Z, p_, T_, sMod_(RInd));

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsNone(CpMCv);

        //- Functions for volumetric coefficient of thermal expansion [1/T]
        scalarFieldFunctionsNone(alphav);

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

#include "linearLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
