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
    Foam::Function1Types::RHWindSpecificDissipationRateProfile

Description
    RHWindSpecificDissipationRateProfile for ABL inlet specification

See also
    Foam::Function1

SourceFiles
    RHWindSpecificDissipationRateProfile.C

\*---------------------------------------------------------------------------*/

#ifndef RHWindSpecificDissipationRateProfile_H
#define RHWindSpecificDissipationRateProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/RHWindProfile/RHWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
           Class RHWindSpecificDissipationRateProfile Declaration
\*---------------------------------------------------------------------------*/

class RHWindSpecificDissipationRateProfile
:
    public RHWindKineticEnergyProfile
{
public:

    // Runtime type information
    TypeName("RHWindSpecificDissipationRateProfile");


    // Constructors

        //- Construct from entry name and dictionary
        RHWindSpecificDissipationRateProfile
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        RHWindSpecificDissipationRateProfile
        (
            const RHWindSpecificDissipationRateProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return
                tmp<Function1<scalar>>
                (
                    new RHWindSpecificDissipationRateProfile(*this)
                );
        }


       //- Destructor
       virtual ~RHWindSpecificDissipationRateProfile();


    // Member Functions

        //- Return value for coordinate z
       virtual scalar value(const scalar z_) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RHWindSpecificDissipationRateProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
