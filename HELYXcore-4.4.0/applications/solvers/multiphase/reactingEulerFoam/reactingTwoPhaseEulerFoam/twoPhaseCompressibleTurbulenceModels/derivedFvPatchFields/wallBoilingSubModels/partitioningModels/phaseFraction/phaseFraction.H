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
    Foam::wallBoilingModels::phaseFraction

Description
    Wall heat-flux partitioned according to the phase volume fraction.

SourceFiles
    phaseFraction.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFraction_H
#define phaseFraction_H

#include "derivedFvPatchFields/wallBoilingSubModels/partitioningModels/partitioningModel/partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class phaseFraction Declaration
\*---------------------------------------------------------------------------*/

class phaseFraction
:
    public partitioningModel
{

public:

    //- Runtime type information
    TypeName("phaseFraction");


    // Constructors

        //- Construct from a dictionary
        phaseFraction(const dictionary& dict);


    //- Destructor
    virtual ~phaseFraction();


    // Member Functions

        //- Calculate and return the wall heat-flux partioning
        virtual tmp<scalarField> fLiquid(const scalarField& alphaLiquid) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
