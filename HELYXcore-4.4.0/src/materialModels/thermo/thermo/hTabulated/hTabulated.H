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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::hTabulated

Description
    Tabulated properties thermodynamics package
    templated into the EquationOfState.

    thermodynamics
    {
        //- Heat capacity table as function of P and T (required)
        CpTableCoeffs
        {
            fileName        "$FOAM_CASE/constant/CpTableName.dat";
            outOfBounds     clamp;
        }

        //- Heat of formation (required)
        Hf 1901;
    }

    The Cp tables are described in terms of
    the pressure and temperature
    For example, the table:

                           temperature (K)
                       |  300   450   600
                  -----+-----------------
                  1e5  |   A     B     C
    pressure (Pa) 1.1e5|   D     E     F
                  2e5  |   G     H     I


    Is specified by the following:

       (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        )

    All pressures and temperatures must be in ascending order
    in each list.

SourceFiles
    hTabulatedI.H
    hTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef hTabulated_H
#define hTabulated_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hTabulated Declaration
\*---------------------------------------------------------------------------*/

class hTabulated
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

        //- 2D look up table Cp = function of pressure and temperature
        //- Specific heat at constant pressure [J/(kg.K)]
        autoPtr<interpolation2DTable<scalar>> CpTable_;

        //- Heat of formation
        scalar hf_;

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
        hTabulated
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<hTabulated> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hTabulated> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hTabulated");

        //- Heat capacity at constant pressure [J/kg/K]
        scalarFieldFunctionsThree(Cp, p_, T_, sMod_(CpContributionInd));

        //- Absolute enthalpy [J/kg]
        scalarFieldFunctionsThree(ha, p_, T_, sMod_(hContributionInd));

        //- Sensible enthalpy [J/kg]
        scalarFieldFunctionsThree(hs, p_, T_, sMod_(hContributionInd));

        //- Entropy [J/kg/K]
        scalarFieldFunctionsThree(s, p_, T_, sMod_(spContributionInd));

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

#include "hTabulatedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
