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
    (c) 2012-2021 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::hPower

Description
    Enthalpy based thermodynamics package using a power function of temperature
    for the constant heat capacity at constant volume which is particularly
    suitable for solids at low temperatures:

    \verbatim
        Cp = c0*pow(T/Tref, n0)
    \endverbatim

Usage
    \table
        Property     | Description
        c0           | Reference heat capacity at constant pressure [J/kg/K]
        n0           | Exponent of the power function
        Tref         | Reference temperature [K]
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hPower thermodynamics for air:

    \verbatim
    thermodynamics hPower;
    thermodynamicsCoeffs
    {
        c0          230;
        Tref        470;
        n0          3;
        Hf          0;
    }
    \endverbatim

SourceFiles
    hPowerI.H
    hPower.C

\*---------------------------------------------------------------------------*/

#ifndef hPower_H
#define hPower_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hPower Declaration
\*---------------------------------------------------------------------------*/

class hPower
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

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar hf_;

        //- Quick model access
        enum
        {
            hsInd,
            hfInd,
            CpContributionInd,
            hContributionInd,
            spContributionInd,
            CpInd,
            CpMCvInd,
            rhoInd,
            haInd,
            modelsEnumSize_
        };


public:


    // Constructors

        //- Construct from dictionary
        hPower
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct by coefficient mixing (for chemistry)
        hPower
        (
            const UPtrList<baseModels<scalar>>& components,
            const scalarField& coeffs,
            const word& name = "linearCombination"
        );

        //- Construct and return a clone
        inline autoPtr<hPower> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hPower> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hPower");

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsTwo(Cp, T_, sMod_(CpContributionInd));

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsTwo(ha, sMod_(hsInd), sMod_(hfInd));

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

#include "hPowerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
