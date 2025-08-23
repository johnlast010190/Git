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
    Foam::muSutherland

Description
    Transport package using Sutherland's formula for viscosity:

    \verbatim
        mu = As*sqrt(T)/(1 + Ts/T)
    \endverbatim


Usage
    \table
        Property     | Description
        As           | Sutherland constant kg/(ms K^1/2)
        Ts           | Sutherland temperature [K]
    \endtable

    Example specification of mu Sutherland for air:

    \verbatim
    muModel muSutherland;
    muModelCoeffs
    {
        As              1.458e-06;
        Ts              110.4;
    }
    \endverbatim

SourceFiles
    muSutherland.C

\*---------------------------------------------------------------------------*/

#ifndef muSutherland_H
#define muSutherland_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class muSutherland Declaration
\*---------------------------------------------------------------------------*/

class muSutherland
:
    public materialModel,
    public muModel
{
    // Private data

        // Sutherland's coefficients
        dimensionedScalar As_;
        dimensionedScalar Ts_;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        dimensionedScalar readCoeff
        (
            const word& coeffName,
            const dimensionSet& dims,
            const dictionary& dict
        );


public:

    // Constructors

        //- Construct from dictionary
        muSutherland
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<muSutherland> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("muSutherland");

        //- Dynamic viscosity [kg/ms]
        scalarFieldFunctionsOne(mu, T_);

        //- Update table
        virtual void updateTable(const word& modelName) override
        {}

        // No dependencies
        materialDepList();

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

#include "muSutherlandI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
