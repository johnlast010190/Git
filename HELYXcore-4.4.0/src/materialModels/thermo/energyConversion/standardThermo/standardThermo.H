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
    (c) 2021-2025 Engys Ltd.

Class
    Foam::standardThermo

SourceFiles
    standardThermo.C

\*---------------------------------------------------------------------------*/

#ifndef standardThermo_H
#define standardThermo_H

#include "materialModels/materialTables/materialTables.H"
#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class standardThermo Declaration
\*---------------------------------------------------------------------------*/

class standardThermo
:
    public materialModel,
    public gammaModel,
    public gModel,
    public aModel
{
    // Private data

        //- Quick model access
        enum {CpMCv, Cp, s, ha, eaInd, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("standardThermo");


    // Constructors

        //- Construct from components
        standardThermo
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- New constructor
        static autoPtr<standardThermo> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    ~standardThermo()
    {}


     // Member functions

        //- Gamma = Cp/Cv []
        scalarFieldFunctionsTwo(gamma, sMod_(Cp), sMod_(CpMCv));

        //- Gibbs free energy [J/kg]
        scalarFieldFunctionsThree(g, T_, sMod_(s), sMod_(ha));

        //- Helmholtz free energy [J/kg]
        scalarFieldFunctionsThree(a, T_, sMod_(s), sMod_(eaInd));

        //- Return model dependence list
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

#include "standardThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
