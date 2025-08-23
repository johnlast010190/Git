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
    (c) 2019-2023 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::rPolynomialLaw

Description
    Reciprocal polynomial equation of state for liquids and solids

    \verbatim
        1/rho = C[0] + C[1]*T + C[2]*sqr(T) - C[3]*p - C[4]*p*T
    \endverbatim

    This polynomial for the reciprocal of the density provides a much better fit
    than the equivalent polynomial for the density and has the advantage that it
    support coefficient mixing to support liquid and solid mixtures in an
    efficient manner.

Usage
    \table
        Property     | Description
        C            | Density polynomial coefficients
    \endtable

    Example specification of the rPolynomial equation of state for pure water:
    \verbatim
    equationOfStateCoeffs
    {
        C (0.001278 -2.1055e-06 3.9689e-09 4.3772e-13 -2.0225e-16);
    }
    \endverbatim
    Note: This fit is based on the small amount of data which is freely
    available for the range 20-65degC and 1-100bar.

SourceFiles
    rPolynomialLawI.H
    rPolynomialLaw.C

\*---------------------------------------------------------------------------*/

#ifndef rPolynomialLaw_H
#define rPolynomialLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class rPolynomialLaw Declaration
\*---------------------------------------------------------------------------*/

class rPolynomialLaw
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
    // Private Classes

        //- Coefficient list class
        class coeffList
        :
            public VectorSpace<coeffList, scalar, 5>
        {
        public:

            // Constructors

                //- Construct null
                inline coeffList()
                {}

                //- Construct from Istream
                inline coeffList(Istream& is)
                :
                    VectorSpace<coeffList, scalar, 5>(is)
                {}
        };


    // Private member data

        //- Quick model access
        enum {R, modelsEnumSize_};

        //- Density coefficients
        coeffList C_;


public:

    // Constructors

        //- Construct from dictionary and object registry
        rPolynomialLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<rPolynomialLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<rPolynomialLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("rPolynomial");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return p_->isConst();
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return T_->isConst() && p_->isConst();
        }

        //- Functions for density [kg/m^3]
        scalarFieldFunctionsTwo(rho, p_, T_);

        //- Functions for enthalpy contribution [J/kg]
        scalarFieldFunctionsTwo(hContribution, p_, T_);

        //- Functions for Cp contribution [J/kg/K]
        scalarFieldFunctionsTwo(CpContribution, p_, T_);

        //- Functions for internal energy contribution [J/kg/K]
        scalarFieldFunctionsTwo(eContribution, p_, T_);

        //- Functions for Cv contribution [J/kg/K]
        scalarFieldFunctionsTwo(CvContribution, p_, T_);

        //- Functions for entropy contribution to the integral of Cp/T [J/kg/K]
        scalarFieldFunctionsNone(spContribution);

        //- Functions for entropy contribution to the integral of Cv/T [J/kg/K]
        scalarFieldFunctionsNone(svContribution);

        //- Functions for compression factor []
        scalarFieldFunctionsThree(Z, p_, T_, sMod_(R));

        //- Functions for compressibility drho/dp [s^2/m^2]
        scalarFieldFunctionsTwo(psi, p_, T_);

        //- Functions for (Cp - Cv) [J/kg/K]
        scalarFieldFunctionsTwo(CpMCv, p_, T_);

        //- Functions for volumetric coefficient of thermal expansion [1/T]
        scalarFieldFunctionsTwo(alphav, p_, T_);

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

#include "rPolynomialLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
