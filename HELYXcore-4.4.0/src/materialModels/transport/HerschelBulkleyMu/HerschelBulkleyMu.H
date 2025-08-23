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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::HerschelBulkleyMu

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    HerschelBulkleyMu.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkleyMu_H
#define HerschelBulkleyMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HerschelBulkleyMu Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkleyMu
:
    public materialModel,
    public muModel
{
    // Private data

        //- Model coefficients
        scalar k_;
        scalar n_;
        scalar tau0_;
        scalar mu0_;

        //- Quick model access
        enum {sRateInd, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("HerschelBulkleyMu");


    // Constructors

        //- Construct from components
        HerschelBulkleyMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~HerschelBulkleyMu()
    {}


    // Member Functions

        //- Dynamic viscosity [kg/ms]
        scalarFieldFunctionsOne(mu, sMod_(sRateInd));

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

#include "HerschelBulkleyMuI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
