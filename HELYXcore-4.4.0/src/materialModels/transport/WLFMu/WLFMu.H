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
    (c) 2020-2021 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::WLFMu

Description
    Transport package using the Williams-Landel-Ferry model for viscosity of
    polymer melts:

    Dynamic viscosity [kg/m.s]
    \verbatim
        mu = mu0*exp(-C1*(T - Tr)/(C2 + T - Tr))
    \endverbatim

    References:
    \verbatim
        Williams, M. L., Landel, R. F., & Ferry, J. D. (1955).
        The temperature dependence of relaxation mechanisms
        in amorphous polymers and other glass-forming liquids.
        Journal of the American Chemical society, 77(14), 3701-3707.
    \endverbatim

    The thermal conductivity is obtained using a constant Prandtl number.

Usage
    \table
        Property     | Description
        mu0          | Reference dynamic viscosity [Pa.s]
        Tref         | Reference temperature [K]
        C1           | WLF constant []
        C2           | WLF constant [K]
        Pr           | Constant Prandtl number []
    \endtable

    Example specification of WLFTransport for a polymer:

    \verbatim
    muModelCoeffs
    {
        mu0         50000;
        Tr          416;
        C1          20.4;
        C2          101.6;
        Pr          10000;
    }
    \endverbatim

SourceFiles
    WLFMuI.H
    WLFMu.C

\*---------------------------------------------------------------------------*/

#ifndef WLFMu_H
#define WLFMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WLFMu Declaration
\*---------------------------------------------------------------------------*/

class WLFMu
:
    public materialModel,
    public muModel
{
    // Private data

        //- Dynamic viscosity at the reference temperature [Pa.s]
        scalar mu0_;

        //- Reference temperature [T]
        scalar Tr_;

        //- WLF coefficient 1 []
        scalar C1_;

        //- WLF coefficient 2 [T]
        scalar C2_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


public:

    //- Runtime type information
    TypeName("WLF");


    // Constructors

        //- Construct from components
        WLFMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~WLFMu()
    {}


    // Member Functions

        //- Dynamic viscosity [kg/ms]
        scalarFieldFunctionsOne(mu, T_);

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

#include "WLFMuI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
