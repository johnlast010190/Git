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
    Foam::hPolynomial

Description
    Enthalpy based thermodynamics package using a polynomial function
    of temperature for the constant heat capacity at constant pressure:

    \verbatim
        Cp = a0 + a1*T + a2*sqr(T) + a3*pow3(T) + a4*pow4(T)
           + a5*pow(T, 5) + a6*pow(T, 6) + a7*pow(T, 7)
    \endverbatim

    The polynomial function is templated on the order of the polynomial which
    defaults to 8.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        CpCoeffs<8>  | Specific heat at constant pressure polynomial coeffs
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics hPolynomial;
    thermodynamicsCoeffs
    {
        Hf              0;
        Sf              0;
        CpCoeffs        8 ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim


Note
    Specific heat at constant pressure polynomial coefficients evaluate to an
    expression in [J/kg/K].

SourceFiles
    hPolynomialI.H
    hPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef hPolynomial_H
#define hPolynomial_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hPolynomial Declaration
\*---------------------------------------------------------------------------*/

class hPolynomial
:
    public coefficientMixture,
    public materialModel,
    public CpModel,
    public haModel,
    public hsModel,
    public hfModel,
    public sModel,
    public dCpdTModel,
    public CvModel,
    public esModel,
    public eaModel
{
    // Private member data

        //- Heat of formation
        scalar hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure polynomial coeffs [J/kg/K/K^i]
        polynomialFunction CpCoeffs_;

        //- enthalpy polynomial coeffs [J/kg/K^i]
        //  Derived from Cp coeffs. Relative to Tstd.
        polynomialFunction hCoeffs_;

        //- Entropy polynomial coeffs [J/kg/K/K^i]
        //  Derived from Cp coeffs. Relative to Tstd.
        polynomialFunction sCoeffs_;

        //- Quick model access
        enum
        {
            CpContributionInd,
            hContributionInd,
            spContributionInd,
            CpInd,
            CpMCvInd,
            rhoInd,
            hsInd,
            haInd,
            modelsEnumSize_
        };


public:


    // Constructors

        //- Construct from dictionary
        hPolynomial
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct by coefficient mixing (for chemistry)
        hPolynomial
        (
            const UPtrList<baseModels<scalar>>& models,
            const scalarField& coeffs,
            const word& name = "linearCombination"
        );

        //- Construct and return a clone
        inline autoPtr<hPolynomial> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hPolynomial> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hPolynomial");

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsTwo(Cp, T_, sMod_(CpContributionInd));

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsTwo(ha, T_, sMod_(hContributionInd));

        //- Sensible enthalpy [J/kg]
        scalarFieldFunctionsTwo(hs, T_, sMod_(hContributionInd));

        //- Entropy [J/kg/K]
        scalarFieldFunctionsTwo(s, T_, sMod_(spContributionInd));

        //- Enthalpy of formation [J/kg]
        scalarFieldFunctionsZero(hf);

        //- Temperature derivative of heat capacity at constant pressure
        scalarFieldFunctionsOne(dCpdT, T_);

        //- Heat capacity at constant volume [J/kg/K]
        scalarFieldFunctionsTwo(Cv, sMod_(CpInd), sMod_(CpMCvInd));

        //- Sensible internal energy [J/kg]
        scalarFieldFunctionsThree(es, p_, sMod_(rhoInd), sMod_(hsInd));

        //- Absolute internal energy [J/kg]
        scalarFieldFunctionsThree(ea, p_, sMod_(rhoInd), sMod_(haInd));

        //- Should be used as premixed gibbs free energy [J/kg]
        virtual scalar gStd(const scalar T) const override;

        //- Delta entropy [J/kg/K]
        virtual scalar deltaS(const scalar T) const override;

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override;

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

#include "hPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
