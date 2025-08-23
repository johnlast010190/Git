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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::kappaSutherland

Description
    Using Sutherland's formula for the thermal conductivity (for gases)
    is obtained from:

    \verbatim
        kappa = mu*Cv*(1.32 + 1.77*R/Cv)
    \endverbatim

Usage
    Example specification of kappa Sutherland for air:

    \verbatim
    kappaModel kappaSutherland;
    \endverbatim

SourceFiles
    kappaSutherland.C

\*---------------------------------------------------------------------------*/

#ifndef kappaSutherland_H
#define kappaSutherland_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class kappaSutherland Declaration
\*---------------------------------------------------------------------------*/

class kappaSutherland
:
    public materialModel,
    public kappaModel
{
    // Private data

        //- Quick model access
        enum {Cv, mu, R, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary
        kappaSutherland
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<kappaSutherland> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("kappaSutherland");

        //- Thermal conductivity [W/mK]
        scalarFieldFunctionsThree(kappa, sMod_(Cv), sMod_(mu), sMod_(R));

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

#include "kappaSutherlandI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
