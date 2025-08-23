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
    (c) 2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::AntoineExtended

Description
    Extended Antoine equation for the vapour pressure.

    \f[
        \log (p) = A + \frac{B}{C + T} + D \log (T) + E T^F
    \f]

    Coefficients \f$A\f$, \f$B\f$, \f$C\f$, \f$D\f$, \f$E\f$ and \f$F\f$ are
    to be supplied and should be suitable for natural logarithms and
    temperatures in Kelvin.

SourceFiles
    AntoineExtended.C

\*---------------------------------------------------------------------------*/

#ifndef AntoineExtended_H
#define AntoineExtended_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AntoineExtended Declaration
\*---------------------------------------------------------------------------*/

class AntoineExtended
:
    public materialModel,
    public pSatModel,
    public pSatPrimeModel,
    public lnPSatModel
{
    // Private data

        //- Constant A
        scalar A_;

        //- Constant B
        scalar B_;

        //- Constant C
        scalar C_;

        //- Constant D
        scalar D_;

        //- Constant F
        scalar F_;

        //- Constant E
        //  (after F so F's dimensions can be used in the construction)
        scalar E_;


public:

    //- Runtime type information
    TypeName("AntoineExtended");


    // Constructors

        //- Construct from components
        AntoineExtended
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~AntoineExtended()
    {}


    // Member Functions

        //- Saturation pressure [Pa]
        scalarFieldFunctionsOne(pSat, T_);

        //- Saturation pressure derivative w.r.t. temperature [Pa/K]
        scalarFieldFunctionsOne(pSatPrime, T_);

        //- Natural log of the saturation pressure
        scalarFieldFunctionsOne(lnPSat, T_);

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

#include "AntoineExtendedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
