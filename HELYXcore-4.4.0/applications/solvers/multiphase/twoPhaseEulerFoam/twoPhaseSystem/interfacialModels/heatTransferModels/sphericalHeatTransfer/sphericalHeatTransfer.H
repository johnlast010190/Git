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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::heatTransferModels::sphericalHeatTransfer

Description
    Model which applies an analytical solution for heat transfer from the
    surface of a sphere to the fluid within the sphere.

SourceFiles
    sphericalHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalHeatTransfer_H
#define sphericalHeatTransfer_H

#include "interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class sphericalHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class sphericalHeatTransfer
:
    public heatTransferModel
{
public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct from components
        sphericalHeatTransfer
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~sphericalHeatTransfer();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
