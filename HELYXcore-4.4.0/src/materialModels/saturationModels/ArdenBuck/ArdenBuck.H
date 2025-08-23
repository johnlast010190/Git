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
    Foam::ArdenBuck

Description
    ArdenBuck equation for the vapour pressure of moist air.

SourceFiles
    ArdenBuck.C

\*---------------------------------------------------------------------------*/

#ifndef ArdenBuck_H
#define ArdenBuck_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ArdenBuck Declaration
\*---------------------------------------------------------------------------*/

class ArdenBuck
:
    public materialModel,
    public pSatModel,
    public pSatPrimeModel,
    public lnPSatModel
{
    static const scalar zeroC;
    static const scalar A;
    static const scalar B;
    static const scalar C;
    static const scalar D;


public:

    //- Runtime type information
    TypeName("ArdenBuck");


    // Constructors

        //- Construct from components
        ArdenBuck
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~ArdenBuck()
    {}


    // Member Functions

        //- Exponent divided by the temperature
        scalar xByTC(scalar TC) const;

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

#include "ArdenBuckI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
