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
    (c) 2025 Engys Ltd.

Class
    Foam::speedSoundGas

Description
    Computes the speed of sound in gases using the ideal gas law:

    \verbatim
        c0 = sqrt(gamma*R*T)
    \endverbatim

Usage
    Example specification of speed of sound model for gases:

    \verbatim
        c0Model speedSoundGas;
    \endverbatim

SourceFiles
    speedSoundGas.C

\*---------------------------------------------------------------------------*/

#ifndef speedSoundGas_H
#define speedSoundGas_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class speedSoundGas Declaration
\*---------------------------------------------------------------------------*/

class speedSoundGas
:
    public materialModel,
    public c0Model
{
    // Private data

        //- Quick model access
        enum {gamma, R, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary
        speedSoundGas
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<speedSoundGas> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("speedSoundGas");

        //- Speed of Sound [m/s]
        scalarFieldFunctionsThree(c0, T_, sMod_(gamma), sMod_(R));

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

#include "speedSoundGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
