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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::matAbsoluteEnthalpy

SourceFiles
    matAbsoluteEnthalpy.C

\*---------------------------------------------------------------------------*/

#ifndef matAbsoluteEnthalpy_H
#define matAbsoluteEnthalpy_H

#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class matAbsoluteEnthalpy Declaration
\*---------------------------------------------------------------------------*/

class matAbsoluteEnthalpy
:
    public materialModel,
    public CpvModel,
    public heModel
{
    // Private data

        //- Quick model access
        enum {Cp, ha, modelsEnumSize_};

public:

    // Constructors

        //- Construct from dictionary
        matAbsoluteEnthalpy
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<matAbsoluteEnthalpy> clone() const;


    //- Destructor
    virtual ~matAbsoluteEnthalpy()
    {};


    // Member functions

        //- Return the instantiated type name
        TypeName("absoluteEnthalpy");

        //- Sensible internal energy type name
        static word name()
        {
            return "ha";
        }

        //- Heat capacity at constant volume [J/kg/K]
        scalarFieldFunctionsOne(Cpv, sMod_(Cp));

        //- Sensible internal energy [J/kg]
        scalarFieldFunctionsOne(he, sMod_(ha));

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
        virtual bool read() override
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "matAbsoluteEnthalpyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
