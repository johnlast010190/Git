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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::kappaPrandtl

Usage
    \table
        Property     | Description
        Pr           | Constant Prandtl number []
    \endtable

    Example specification:
    \verbatim
    kappaModelCoeffs
    {
        Pr              0.71;
    }
    \endverbatim

SourceFiles
    kappaPrandtl.C

\*---------------------------------------------------------------------------*/

#ifndef kappaPrandtl_H
#define kappaPrandtl_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class kappaPrandtl Declaration
\*---------------------------------------------------------------------------*/

class kappaPrandtl
:
    public materialModel,
    public kappaModel
{
    // Private data

        //- Reciprocal Prandtl Number []
        scalar rPr_;

        //- Quick model access
        enum {Cp, mu, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary
        kappaPrandtl
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<kappaPrandtl> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("kappaPrandtl");

        //- Is it isotropic model?
        virtual bool isotropic() const override
        {
            return true;
        }

        //- Thermal conductivity [W/mK]
        scalarFieldFunctionsTwo(kappa, sMod_(Cp), sMod_(mu));

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

#include "kappaPrandtlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //