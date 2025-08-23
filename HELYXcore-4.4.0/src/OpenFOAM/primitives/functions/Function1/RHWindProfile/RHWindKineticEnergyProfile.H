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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::Function1Types::RHWindKineticEnergyProfile

Description
    RHWindKineticEnergyProfile (RichardsHoxey) for ABL inlet specification
    constant with z
See also
    Foam::Function1

SourceFiles
    RHWindKineticEnergyProfile.C

\*---------------------------------------------------------------------------*/

#ifndef RHWindKineticEnergyProfile_H
#define RHWindKineticEnergyProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/RHWindProfile/RHWindProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                Class RHWindKineticEnergyProfile Declaration
\*---------------------------------------------------------------------------*/

class RHWindKineticEnergyProfile
:
    public RHWindProfile
{
public:

    // Runtime type information
    TypeName("RHWindKineticEnergyProfile");


    // Constructors

        //- Construct from entry name and dictionary
        RHWindKineticEnergyProfile
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        RHWindKineticEnergyProfile
        (
            const RHWindKineticEnergyProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return
                tmp<Function1<scalar>>
                (
                    new RHWindKineticEnergyProfile(*this)
                );
        }

       //- Destructor
       virtual ~RHWindKineticEnergyProfile();


    // Member Functions

        //- Return value for coordinate z
       virtual scalar value(const scalar z_) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RHWindKineticEnergyProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
