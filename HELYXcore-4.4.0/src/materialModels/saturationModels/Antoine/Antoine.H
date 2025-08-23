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
    Foam::Antoine

Description
    Antoine equation for the vapour pressure.

    \f[
        \log p = A + \frac{B}{C + T}
    \f]

    Coefficients \f$A\f$, \f$B\f$ and \f$C\f$ are to be supplied and should be
    suitable for natural logarithms and temperatures in Kelvin.

SourceFiles
    Antoine.C

\*---------------------------------------------------------------------------*/

#ifndef Antoine_H
#define Antoine_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Antoine Declaration
\*---------------------------------------------------------------------------*/

class Antoine
:
    public materialModel,
    public pSatModel,
    public pSatPrimeModel,
    public lnPSatModel,
    public TsatModel
{
    // Private data

        //- Constant A
        scalar A_;

        //- Constant B
        scalar B_;

        //- Constant C
        scalar C_;


public:

    //- Runtime type information
    TypeName("Antoine");


    // Constructors

        //- Construct from components
        Antoine
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~Antoine()
    {}


    // Member Functions

        //- Saturation pressure [Pa]
        scalarFieldFunctionsOne(pSat, T_);

        //- Saturation pressure derivative w.r.t. temperature [Pa/K]
        scalarFieldFunctionsOne(pSatPrime, T_);

        //- Natural log of the saturation pressure
        scalarFieldFunctionsOne(lnPSat, T_);

        //- Saturation temperature [K]
        scalarFieldFunctionsOne(Tsat, p_);

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

#include "AntoineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
