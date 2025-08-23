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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::wallBoilingModels::Lavieville

Description
    Lavieville wall heat flux partitioning model.

    Model parameters:
        alphaCrit: critical liquid fraction

    Reference:
    \verbatim
        Lavieville, J., Quemerais, E., Mimouni, S., Boucker, M., &
        Mechitoua, N. (2006).
        NEPTUNE CFD V1. 0 theory manual.
        NEPTUNE report Nept_2004_L1, 2(3).
    \endverbatim

SourceFiles
    Lavieville.C

\*---------------------------------------------------------------------------*/

#ifndef Lavieville_H
#define Lavieville_H

#include "derivedFvPatchFields/wallBoilingSubModels/partitioningModels/partitioningModel/partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class Lavieville Declaration
\*---------------------------------------------------------------------------*/

class Lavieville
:
    public partitioningModel
{
    // Private data

        //- Critical liquid fraction
        scalar alphaCrit_;


public:

    //- Runtime type information
    TypeName("Lavieville");


    // Constructors

        //- Construct from a dictionary
        Lavieville(const dictionary& dict);


    //- Destructor
    virtual ~Lavieville();


    // Member Functions

        //- Liquid blending function
        virtual tmp<scalarField> fLiquid(const scalarField& alphaLiquid) const;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
