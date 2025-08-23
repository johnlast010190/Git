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
    Foam::Function1Types::AIJWindProfile

Description
    AIJWindProfile for ABL inlet specification
    Usage:
    \verbatim
        <entryName> <AIJWindProfileFunction>;
    <entryName>Coeffs
    {
        Cmu 0.09;
        Uref 5;
        alpha 0.2;
        Href 10;
        zG 10;
    }
    \endverbatim

See also
    Foam::Function1

SourceFiles
    AIJWindProfile.C

\*---------------------------------------------------------------------------*/

#ifndef AIJWindProfile_H
#define AIJWindProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class AIJWindProfile Declaration
\*---------------------------------------------------------------------------*/

class AIJWindProfile
:
    public Function1<scalar>
{
protected:

    // Protected Member Data

        //- Cmu constant
        scalar Cmu_;

        //- Reference velocity
        scalar Uref_;

        //- Reference height
        scalar Href_;

        //- Power law exponent determined by terrain category
        scalar alpha_;

        //- boundary layer height determined by terrain category
        scalar zG_;

    // Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    // Runtime type information
    TypeName("AIJWindProfile");


    // Constructors

        //- Construct from entry name and dictionary
        AIJWindProfile(const word& entryName, const dictionary& dict);

        //- Copy constructor
        AIJWindProfile(const AIJWindProfile& wp);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return tmp<Function1<scalar>>(new AIJWindProfile(*this));
        }


    //- Destructor
    virtual ~AIJWindProfile();


    // Member Functions

        //- Return value for coordinate z
        virtual scalar value(const scalar z) const override;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const AIJWindProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
