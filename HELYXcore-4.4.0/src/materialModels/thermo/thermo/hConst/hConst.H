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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::hConst

Description
    Enthalpy based thermodynamics package using a constant heat capacity
    at constant pressure.

    The sensible enthalpy is evaluated as:
    \verbatim
        hs = Cp*(T - Tref) + Hsref
    \endverbatim

Usage
    \table
        Property     | Description
        Cp           | Constant Heat capacity at constant pressure [J/kg/K]
        Tref         | Reference temperature [K] (defaults to 0)
        Hsref        | Reference sensible enthalpy [J/kg] (defaults to 0)
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of hConst for air:

    \verbatim
    thermodynamics hConst;
    thermodynamicsCoeffs
    {
        Cp          1007;
        Hf          0;
    }
    \endverbatim

SourceFiles
    hConstI.H
    hConst.C

\*---------------------------------------------------------------------------*/

#ifndef hConst_H
#define hConst_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "mixtures/coefficientMixture/coefficientMixture.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class hConst Declaration
\*---------------------------------------------------------------------------*/

class hConst
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

        //- Heat capacity [J/kg/K]
        scalar Cp_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Reference temperature around which to linearise [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar hsref_;

        //- Entropy [J/kg/K]
        scalar Sf_;

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
        hConst
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct by coefficient mixing (for chemistry)
        hConst
        (
            const UPtrList<baseModels<scalar>>& models,
            const scalarField& coeffs,
            const word& name = "linearCombination"
        );

        //- Construct and return a clone
        inline autoPtr<hConst> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hConst> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hConst");

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsOne(Cp, sMod_(CpContributionInd));

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

#include "hConstI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
