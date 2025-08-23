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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::BirdCarreauMu

Description
    An incompressible Bird-Carreau non-Newtonian viscosity model.

    The Bird-Carreau-Yasuda form is also supported if the optional "a"
    coefficient is specified.  "a" defaults to 2 for the Bird-Carreau model.

SourceFiles
    BirdCarreauMu.C

\*---------------------------------------------------------------------------*/

#ifndef BirdCarreauMu_H
#define BirdCarreauMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BirdCarreauMu Declaration
\*---------------------------------------------------------------------------*/

class BirdCarreauMu
:
    public materialModel,
    public muModel
{
    // Private data

        //- Bird-carreau coefficients
        scalar mu0_;
        scalar muInf_;
        scalar k_;
        scalar n_;
        scalar a_;

        //- Quick model access
        enum {sRateInd, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("BirdCarreauMu");


    // Constructors

        //- Construct from dictionary
        BirdCarreauMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~BirdCarreauMu()
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

#include "BirdCarreauMuI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

