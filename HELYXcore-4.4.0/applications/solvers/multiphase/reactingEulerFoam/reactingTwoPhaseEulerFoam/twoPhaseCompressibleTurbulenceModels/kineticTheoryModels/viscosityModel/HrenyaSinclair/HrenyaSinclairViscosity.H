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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::kineticTheoryModels::viscosityModels::HrenyaSinclair

Description

SourceFiles
    HrenyaSinclairViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef HrenyaSinclair_H
#define HrenyaSinclair_H

#include "kineticTheoryModels/viscosityModel/viscosityModel/viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HrenyaSinclair Declaration
\*---------------------------------------------------------------------------*/

class HrenyaSinclair
:
    public viscosityModel
{
    // Private data

        dictionary coeffDict_;

        //- Characteristic length of geometry
        dimensionedScalar L_;


public:

    //- Runtime type information
    TypeName("HrenyaSinclair");


    // Constructors

        //- Construct from components
        HrenyaSinclair(const dictionary& dict);


    //- Destructor
    virtual ~HrenyaSinclair();


    // Member functions

        tmp<volScalarField> nu
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const dimensionedScalar& e
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
