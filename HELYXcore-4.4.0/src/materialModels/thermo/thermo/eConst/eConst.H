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
    Foam::eConst

Description
    Internal energy based thermodynamics package using a constant heat capacity
    at constant volume.

    The sensible internal energy is evaluated as:
    \verbatim
        es = Cv*(T - Tref) + Esref
    \endverbatim

Usage
    \table
        Property     | Description
        Cv           | Constant Heat capacity at constant volume [J/kg/K]
        Tref         | Reference temperature [K] (defaults to 0)
        Esref        | Reference sensible internal energy [J/kg] (defaults to 0)
        Hf           | Heat of formation [J/kg]
    \endtable

    Example specification of eConst for air:

    \verbatim
    thermodynamics eConst;
    thermodynamicsCoeffs
    {
        Cv          724;
        Hf          0;
    }
    \endverbatim


SourceFiles
    eConstI.H
    eConst.C

\*---------------------------------------------------------------------------*/

#ifndef eConst_H
#define eConst_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class eConst Declaration
\*---------------------------------------------------------------------------*/

class eConst
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

        //- Heat capacity [J/kg/K]
        scalar Cv_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Reference temperature [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar esref_;

        //- Quick model access
        enum
        {
            CpMCvInd,
            CpInd,
            CvInd,
            CvContributionInd,
            eContributionInd,
            svContributionInd,
            rhoInd,
            esInd,
            eaInd,
            modelsEnumSize_
        };


public:


    // Constructors

        //- Construct from dictionary
        eConst
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<eConst> clone() const;

        //- Selector from dictionary
        inline static autoPtr<eConst> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("eConst");

        //- Heat capacity at constant volume [J/kg/K]
        scalarFieldFunctionsOne(Cv, sMod_(CvContributionInd));

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
        scalarFieldFunctionsThree(s, T_, sMod_(svContributionInd), sMod_(CpInd));

        //- Enthalpy of formation [J/kg]
        scalarFieldFunctionsZero(hf);

        //- Temperature derivative of heat capacity at constant pressure
        scalarFieldFunctionsOne(dCpdT, T_);

        //- Gibbs free energy of the mixture in the standard state [J/kg]
        // virtual scalar gStd(const scalar T) const override;

        // //- Delta entropy [J/kg/K]
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

#include "eConstI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
