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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::strainRateLaw

Description
    Strain rate computed from gradient: mag(symm(grad(U))). Gradient is cashed
    in most cases. The class, however has it's own handeling of grad(U) update
    in order to make sure it doesn't happen every single call.

SourceFiles
    strainRateLaw.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateLaw_H
#define strainRateLaw_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class strainRateLaw Declaration
\*---------------------------------------------------------------------------*/

class strainRateLaw
:
    public materialModel,
    public strainRateModel
{
    // Private data

        //- Velocity
        // Note: This needs to be pointer to U directly in order to check for
        // updates (cashed gradU).
        const volVectorField* U_;

        //- Gradient grad(U)
        mutable autoPtr<volTensorField> gradU_;

        //- Updating gradient U
        const volTensorField& gradU() const;


public:

    // Constructors

        //- Construct from dictionary
        strainRateLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<strainRateLaw> clone() const;


    // Destructor
    ~strainRateLaw()
    {
        // Avoid destruction if grad(U) is owned by object registry
        if (gradU_.valid() && gradU_().ownedByRegistry())
        {
            gradU_.release();
        }
    }

    // Member functions

        //- Return the instantiated type name
        TypeName("strainRate");

        bool isstrainRateConst() const override
        {
            return false;
        }

        //- Strai rate [1/s]
        declareScalarFieldFunctions(strainRate);

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

        //- Update fields
        virtual bool updateVectorField
        (
            const word& fieldName,
            const volVectorField& volField
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
