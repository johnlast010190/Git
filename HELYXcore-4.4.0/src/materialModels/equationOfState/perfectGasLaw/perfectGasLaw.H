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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::perfectGasLaw

Description
    Perfect gas equation of state:

    \verbatim
        rho = p/(R*T)
    \endverbatim

Usage
    The gas constant R used by the perfect gas equation of state is obtained
    directly from the specie or mixture molecular weight so there is no need to
    provide an equation of \c equationOfState enry in the \c mixture
    specification.

    Example specification of the linear equation of state:

    \verbatim
    equationOfState perfectGas;
    \endverbatim

SourceFiles
    perfectGasLawI.H
    perfectGasLaw.C

\*---------------------------------------------------------------------------*/

#ifndef perfectGasLaw_H
#define perfectGasLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class perfectGasLaw Declaration
\*---------------------------------------------------------------------------*/

class perfectGasLaw
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

        //- Quick model access
        enum {RInd, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary and object registry
        perfectGasLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<perfectGasLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<perfectGasLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("perfectGas");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override;

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsThree(rho, p_, T_, sMod_(RInd));

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsNone(hContribution);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsNone(CpContribution);

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsNone(eContribution);

        //- Functions for Cv contribution [J/kg/K]
        scalarFieldFunctionsNone(CvContribution);

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsTwo(spContribution, p_, sMod_(RInd));

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsNone(svContribution);

        //- Functions for compression factor []
        scalarFieldFunctionsZero(Z);

        //- Functions for compressibility drho/dp [s^2/m^2]
        scalarFieldFunctionsTwo(psi, T_, sMod_(RInd));

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsOne(CpMCv, sMod_(RInd));

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

#include "perfectGasLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
