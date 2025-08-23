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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019-2025 Engys Ltd.

Class
    Foam::hFunctionT

Description
    Thermodynamics package templated on the equation of state, using Function1
    functions for \c cp, \c h and \c s.

    Functions for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        Cp           | Specific heat at constant pressure function
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        Cp              polynomial ( ( 1000 0 ) (-0.05 1) (0.003 2) );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            Cp = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/(kg.K)].

SourceFiles
    hFunctionTI.H
    hFunctionT.C

\*---------------------------------------------------------------------------*/

#ifndef hFunctionT_H
#define hFunctionT_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
#include "primitives/functions/Function1/Constant/Constant.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hFunctionT Declaration
\*---------------------------------------------------------------------------*/

class hFunctionT
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

        //- Specific heat at constant pressure
        autoPtr<Function1<scalar>> Cp_;

        //- Specific heat at constant pressure
        autoPtr<Function1<scalar>> h_;

        //- Quick model access
        enum
        {
            haInd,
            hfInd,
            CpContributionInd,
            hContributionInd,
            spContributionInd,
            CpInd,
            CpMCvInd,
            rhoInd,
            hsInd,
            modelsEnumSize_
        };


public:


    // Constructors

        //- Construct from dictionary
        hFunctionT
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<hFunctionT> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hFunctionT> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hFunctionT");

        virtual bool isCpConst() const override
        {
            return
                (
                    Cp_().type() == Function1Types::Constant<scalar>::typeName
                 || T_->isConst()
                )
             && isModelConst<CpContributionModel>();
        }

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsNoConstTwo(Cp, T_, sMod_(CpContributionInd));

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

#include "hFunctionTI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
