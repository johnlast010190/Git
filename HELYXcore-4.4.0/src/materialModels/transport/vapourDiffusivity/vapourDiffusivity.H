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
    (c) 2025 Engys Ltd.

Class
    Foam::vapourDiffusivity

Description
    Vapour diffusivity for liquid properties.

SourceFiles
    vapourDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef vapourDiffusivity_H
#define vapourDiffusivity_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Function2/APIdiffCoef/APIdiffCoef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class vapourDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class vapourDiffusivity
:
    public materialModel,
    public DModel
{
    // Private data

        //- Vapour diffusivity [m2/s]
        autoPtr<Function2s::APIdiffCoef> D_;

        //- Vapour diffusivity [m2/s]
        autoPtr<Function1<scalar>> Dfunc1_;


public:

    // Constructors

        //- Construct from dictionary
        vapourDiffusivity
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<vapourDiffusivity> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("vapourDiffusivity");

        //- Thermal conductivity [W/mK]
        scalarFieldFunctionsTwo(D, p_, T_);

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override
        {}

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

#include "vapourDiffusivityI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //