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
    Foam::janaf

Description
    Enthalpy based thermodynamics package using JANAF tables:

    \verbatim
        Cp = (((a4*T + a3)*T + a2)*T + a1)*T + a0
        ha = ((((a4/5*T + a3/4)*T + a2/3)*T + a1/2)*T + a0)*T + a5
    \endverbatim

Usage
    \table
        Property     | Description
        Tlow         | Lower temperature limit [K]
        Thigh        | Upper temperature limit [K]
        Tcommon      | Transition temperature from low to high polynomials [K]
        lowCpCoeffs  | Low temperature range heat capacity coefficients
        highCpCoeffs | High temperature range heat capacity coefficients
    \endtable

    Example specification of janafThermo for air:

    \verbatim
    thermodynamics janaf;
    thermodynamicsCoeffs
    {
        Tlow            100;
        Thigh           10000;
        Tcommon         1000;
        lowCpCoeffs
        (
            3.5309628
            -0.0001236595
            -5.0299339e-07
            2.4352768e-09
            -1.4087954e-12
            -1046.9637
            2.9674391
        );
        highCpCoeffs
        (
            2.9525407
            0.0013968838
            -4.9262577e-07
            7.8600091e-11
            -4.6074978e-15
            -923.93753
            5.8718221
        );
    }
    \endverbatim

SourceFiles
    janafI.H
    janaf.C

\*---------------------------------------------------------------------------*/

#ifndef janaf_H
#define janaf_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class janaf Declaration
\*---------------------------------------------------------------------------*/

class janaf
:
    public coefficientMixture,
    public materialModel,
    public CpModel,
    public haModel,
    public hsModel,
    public hfModel,
    public sModel,
    public dCpdTModel,
    public limitModel,
    public CvModel,
    public esModel,
    public eaModel
{
public:

    // Public data

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private member data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;

        //- R
        scalar R_;

        bool isConverted_;

        //- Quick model access
        enum
        {
            CpContribution,
            hContributionInd,
            spContributionInd,
            haInd,
            hfInd,
            RInd,
            CpInd,
            CpMCvInd,
            rhoInd,
            hsInd,
            modelsEnumSize_
        };


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;

public:


    // Constructors

        //- Construct from dictionary
        janaf
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct by coefficient mixing (for chemistry)
        janaf
        (
            const UPtrList<baseModels<scalar>>& models,
            const scalarField& coeffs,
            const word& name = "linearCombination"
        );

        //- Construct and return a clone
        inline autoPtr<janaf> clone() const;

        //- Selector from dictionary
        inline static autoPtr<janaf> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("janaf");

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        scalarFieldFunctionsOne(limit, T_);

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsTwo(Cp, T_, sMod_(CpContribution));

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsTwo(ha, T_, sMod_(hContributionInd));

        //- Entropy [J/kg/K]
        scalarFieldFunctionsTwo(s, T_, sMod_(spContributionInd));

        //- Sensible enthalpy [J/kg]
        scalarFieldFunctionsTwo(hs, sMod_(haInd), sMod_(hfInd));

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

#include "janafI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
