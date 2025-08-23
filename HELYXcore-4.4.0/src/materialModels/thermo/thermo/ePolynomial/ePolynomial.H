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
    (c) 2020 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Description
    Internal energy based thermodynamics package using a polynomial function
    of temperature for the constant heat capacity at constant volume:

    \verbatim
        Cv = CvCoeffs[0] + CvCoeffs[1]*T + CvCoeffs[2]*sqr(T)
           + CvCoeffs[3]*pow3(T) + CvCoeffs[4]*pow4(T)
           + CvCoeffs[5]*pow(T, 5) + CvCoeffs[6]*pow(T, 6)
           + CvCoeffs[7]*pow(T, 7)
    \endverbatim

    The polynomial function is templated on the order of the polynomial which
    defaults to 8.

Usage
    \table
        Property     | Description
        Hf           | Heat of formation
        sf           | Standard entropy
        CvCoeffs<8>  | Specific heat at constant volume polynomial coeffs
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        sf              0;
        CvCoeffs       (1000 -0.05 0.003 0 0 0 0 0);
    }
    \endverbatim

    Note:
        Specific heat at constant volume polynomial coefficients evaluate to an
        expression in [J/kg/K].

SourceFiles
    ePolynomialThermoI.H
    ePolynomialThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef ePolynomial_H
#define ePolynomial_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ePolynomial Declaration
\*---------------------------------------------------------------------------*/

class ePolynomial
:
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

        //- Specific heat at constant volume polynomial coeffs [J/kg/K/K^i]
        polynomialFunction CvCoeffs_;

        //- Internal energy polynomial coeffs [J/kg/K^i]
        //  Derived from Cv coeffs. Relative to Tstd.
        polynomialFunction eCoeffs_;

        //- Entropy polynomial coeffs [J/kg/K/K^i]
        //  Derived from Cv coeffs. Relative to Tstd.
        polynomialFunction sCoeffs_;

        //- Quick model access
        enum
        {
            CvContributionInd,
            eContributionInd,
            svContributionInd,
            CvInd,
            CpInd,
            CpMCvInd,
            rhoInd,
            esInd,
            eaInd,
            modelsEnumSize_
        };


public:


    // Constructors

        //- Construct from dictionary
        ePolynomial
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct by coefficient mixing (for chemistry)
        ePolynomial
        (
            const UPtrList<baseModels<scalar>>& models,
            const scalarField& coeffs,
            const word& name = "linearCombination"
        );

        //- Construct and return a clone
        inline autoPtr<ePolynomial> clone() const;

        //- Selector from dictionary
        inline static autoPtr<ePolynomial> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("ePolynomial");

        //- Heat capacity at constant volume [J/kg/K]
        scalarFieldFunctionsTwo(Cv, T_, sMod_(CvContributionInd));

        //- Sensible internal energy [J/kg]
        scalarFieldFunctionsTwo(es, T_, sMod_(eContributionInd));

        //- Absolute internal energy [J/kg]
        scalarFieldFunctionsTwo(ea, T_, sMod_(eContributionInd));

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsTwo(Cp, sMod_(CvInd), sMod_(CpMCvInd));

        //- Sensible enthalpy [J/kg]
        scalarFieldFunctionsThree(hs, p_, sMod_(esInd), sMod_(rhoInd));

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsThree(ha, p_, sMod_(eaInd), sMod_(rhoInd));

        //- Entropy [J/kg/K]
        scalarFieldFunctionsTwo(s, T_, sMod_(svContributionInd));

        //- Enthalpy of formation [J/kg]
        scalarFieldFunctionsZero(hf);

        //- Temperature derivative of heat capacity at constant pressure
        scalarFieldFunctionsOne(dCpdT, T_);

        //- Should be used as premixed gibbs free energy [J/kg]
        // virtual scalar gStd(const scalar T) const override;

        //- Delta entropy [J/kg/K]
        // virtual scalar deltaS(const scalar T) const override;

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

#include "ePolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
